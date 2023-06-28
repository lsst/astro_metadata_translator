# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

"""Represent a collection of translated headers."""

from __future__ import annotations

__all__ = ("ObservationGroup",)

import logging
from collections.abc import Callable, Iterable, Iterator, MutableMapping, MutableSequence
from typing import TYPE_CHECKING, Any

from .observationInfo import ObservationInfo

if TYPE_CHECKING:
    from .translator import MetadataTranslator

log = logging.getLogger(__name__)


class ObservationGroup(MutableSequence):
    """A collection of `ObservationInfo` headers.

    Parameters
    ----------
    members : iterable of `ObservationInfo` or `dict`-like
        `ObservationInfo` to seed the group membership.  If `dict`-like
        values are used they will be passed to the `ObservationInfo`
        constructor.
    translator_class : `MetadataTranslator`-class, optional
        If any of the members is not an `ObservationInfo`, translator class
        to pass to the `ObservationInfo` constructor.  If `None` the
        translation class will be determined automatically.
    pedantic : `bool`, optional
        If any of the members is not an `ObservationInfo`, passed to the
        `ObservationInfo` constructor to control whether
        a failed translation is fatal or not. `None` indicates that the
        `ObservationInfo` constructor default should be used.
    """

    def __init__(
        self,
        members: Iterable[ObservationInfo | MutableMapping[str, Any]],
        translator_class: type[MetadataTranslator] | None = None,
        pedantic: bool | None = None,
    ) -> None:
        self._members = [
            self._coerce_value(m, translator_class=translator_class, pedantic=pedantic) for m in members
        ]

        # Cache of members in time order
        self._sorted: list[ObservationInfo] | None = None

    def __len__(self) -> int:
        return len(self._members)

    def __delitem__(self, index: int) -> None:  # type: ignore
        del self._members[index]
        self._sorted = None

    def __getitem__(self, index: int) -> ObservationInfo:  # type: ignore
        return self._members[index]

    def __str__(self) -> str:
        results = []
        for obs_info in self._members:
            results.append(f"({obs_info.instrument}, {obs_info.datetime_begin})")
        return "[" + ", ".join(results) + "]"

    def _coerce_value(
        self,
        value: ObservationInfo | MutableMapping[str, Any],
        translator_class: type[MetadataTranslator] | None = None,
        pedantic: bool | None = None,
    ) -> ObservationInfo:
        """Given a value, ensure it is an `ObservationInfo`.

        Parameters
        ----------
        value : `ObservationInfo` or `dict`-like
            Either an `ObservationInfo` or something that can be passed to
            an `ObservationInfo` constructor.
        translator_class : `MetadataTranslator`-class, optional
            If value is not an `ObservationInfo`, translator class to pass to
            the `ObservationInfo` constructor.  If `None` the
            translation class will be determined automatically.
        pedantic : `bool`, optional
            If value is not an `ObservationInfo`, passed to the
            `ObservationInfo` constructor to control whether
            a failed translation is fatal or not. `None` indicates that the
            `ObservationInfo` constructor default should be used.

        Raises
        ------
        ValueError
            Raised if supplied value is not an `ObservationInfo` and can
            not be turned into one.
        """
        if value is None:
            raise ValueError("An ObservationGroup cannot contain 'None'")

        if not isinstance(value, ObservationInfo):
            try:
                kwargs: dict[str, Any] = {"translator_class": translator_class}
                if pedantic is not None:
                    kwargs["pedantic"] = pedantic
                value = ObservationInfo(value, **kwargs)
            except Exception as e:
                raise ValueError("Could not convert value to ObservationInfo") from e

        return value

    def __iter__(self) -> Iterator[ObservationInfo]:
        return iter(self._members)

    def __eq__(self, other: Any) -> bool:
        """Check equality with another group.

        Compare equal if all the members are equal in the same order.
        """
        if not isinstance(other, ObservationGroup):
            return NotImplemented

        for info1, info2 in zip(self, other):
            if info1 != info2:
                return False
        return True

    def __setitem__(  # type: ignore
        self, index: int, value: ObservationInfo | MutableMapping[str, Any]
    ) -> None:
        """Store item in group.

        Item must be an `ObservationInfo` or something that can be passed
        to an `ObservationInfo` constructor.
        """
        value = self._coerce_value(value)
        self._members[index] = value
        self._sorted = None

    def insert(self, index: int, value: ObservationInfo | MutableMapping[str, Any]) -> None:
        value = self._coerce_value(value)
        self._members.insert(index, value)
        self._sorted = None

    def reverse(self) -> None:
        self._members.reverse()

    def sort(self, key: Callable | None = None, reverse: bool = False) -> None:
        self._members.sort(key=key, reverse=reverse)
        if key is None and not reverse and self._sorted is None:
            # Store sorted order in cache
            self._sorted = self._members.copy()

    def extremes(self) -> tuple[ObservationInfo, ObservationInfo]:
        """Return the oldest observation in the group and the newest.

        If there is only one member of the group, the newest and oldest
        can be the same observation.

        Returns
        -------
        oldest : `ObservationInfo`
            Oldest observation.
        newest : `ObservationInfo`
            Newest observation.
        """
        if self._sorted is None:
            self._sorted = sorted(self._members)
        return self._sorted[0], self._sorted[-1]

    def newest(self) -> ObservationInfo:
        """Return the newest observation in the group.

        Returns
        -------
        newest : `ObservationInfo`
            The newest `ObservationInfo` in the `ObservationGroup`.
        """
        return self.extremes()[1]

    def oldest(self) -> ObservationInfo:
        """Return the oldest observation in the group.

        Returns
        -------
        oldest : `ObservationInfo`
            The oldest `ObservationInfo` in the `ObservationGroup`.
        """
        return self.extremes()[0]

    def property_values(self, property: str) -> set[Any]:
        """Return a set of values associated with the specified property.

        Parameters
        ----------
        property : `str`
            Property of an `ObservationInfo`

        Returns
        -------
        values : `set`
            All the distinct values for that property within this group.
        """
        return {getattr(obs_info, property) for obs_info in self}

    def to_simple(self) -> list[MutableMapping[str, Any]]:
        """Convert the group to simplified form.

        Returns
        -------
        simple : `list` of `dict`
            Simple form is a list containing the simplified dict form of
            each `ObservationInfo`.
        """
        return [obsinfo.to_simple() for obsinfo in self]

    @classmethod
    def from_simple(cls, simple: list[dict[str, Any]]) -> ObservationGroup:
        """Convert simplified form back to `ObservationGroup`.

        Parameters
        ----------
        simple : `list` of `dict`
            Object returned by `to_simple`.

        Returns
        -------
        group : `ObservationGroup`
            Reconstructed group.
        """
        return cls(ObservationInfo.from_simple(o) for o in simple)
