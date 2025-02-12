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

from __future__ import annotations

__all__ = ("main",)

import importlib
import logging
import os
from collections.abc import Sequence
from importlib.metadata import entry_points

import click

from ..bin.translate import translate_or_dump_headers
from ..bin.writeindex import write_index_files
from ..bin.writesidecar import write_sidecar_files

# Default regex for finding data files
re_default = r"\.fit[s]?\b"

log = logging.getLogger("astro_metadata_translator")

PACKAGES_VAR = "METADATA_TRANSLATORS"

hdrnum_option = click.option(
    "-n",
    "--hdrnum",
    default=-1,
    help="HDU number to read. If the HDU can not be found, a warning is issued but"
    " reading is attempted using the primary header. The primary header is"
    " always read and merged with this header. Negative number (the default) "
    " indicates that the second header will be merged if the FITS file supports"
    " extended FITS.",
)
regex_option = click.option(
    "-r",
    "--regex",
    default=re_default,
    help="When looking in a directory, regular expression to use to determine whether"
    f" a file should be examined. Default: '{re_default}'",
)
content_option = click.option(
    "-c",
    "--content",
    default="translated",
    type=click.Choice(["translated", "metadata"], case_sensitive=False),
    help="Content to store in JSON file. Options are: "
    "'translated' stores translated metadata in the file; "
    "'metadata' stores raw FITS headers in the file.",
)


@click.group(
    name="astrometadata",
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)
@click.option(
    "--log-level",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False),
    default="INFO",
    help="Python logging level to use.",
)
@click.option(
    "--traceback/--no-traceback", default=False, help="Give detailed trace back when any errors encountered."
)
@click.option(
    "-p",
    "--packages",
    multiple=True,
    help="Python packages or plugin names to import to register additional translators. This is in addition"
    f" to any packages specified in the {PACKAGES_VAR} environment variable (colon-separated"
    " python module names or plugin names).",
)
@click.option(
    "--list-plugins/--no-list-plugins",
    default=False,
    help="List all available registered plugins. If true, the command will return immediately.",
)
@click.pass_context
def main(
    ctx: click.Context, log_level: int, traceback: bool, packages: Sequence[str], list_plugins: bool
) -> None:
    """Execute main click command-line."""
    ctx.ensure_object(dict)

    # Currently we set the log level globally for all Python loggers
    # rather than having a metadata translator logger that all metadata
    # translators can use as a parent. This can potentially cause spurious
    # messages from numexpr. Try to hide those by setting the numexpr env var.
    os.environ["NUMEXPR_MAX_THREADS"] = "8"

    logging.basicConfig(level=log_level)

    # Traceback needs to be known to subcommands
    ctx.obj["TRACEBACK"] = traceback

    plugins = {p.name: p for p in entry_points(group="astro_metadata_translators")}
    if list_plugins:
        from astro_metadata_translator import MetadataTranslator

        print("Builtin translators:")
        for t in sorted(MetadataTranslator.translators):
            print(f"- {t}")
        if plugins:
            print("Available translator plugins grouped by label (use '-p <label>' to activate):")
        for label in sorted(plugins):
            print(f"* {label}:")
            try:
                func = plugins[label].load()
            except Exception as e:
                print(f"  - Unable to load plugin [{e}]")
                continue
            translators = func()
            for t in translators:
                print(f"  - {t}")
        # Exit early with good status.
        raise click.exceptions.Exit(0)

    if ctx.invoked_subcommand is None:
        # Print the help if we were invoked without a subcommand.
        click.echo(ctx.get_help())
        raise click.exceptions.Exit(0)

    packages_set = set(packages)
    if PACKAGES_VAR in os.environ:
        new_packages = os.environ[PACKAGES_VAR].split(":")
        packages_set.update(new_packages)

    # Process import requests
    for m in packages_set:
        if m in plugins:
            try:
                # Loading is sufficient to register the translator.
                plugins[m].load()
            except Exception as e:
                log.warning("Failed to import plugin %s: %s", m, e)
            continue
        try:
            importlib.import_module(m)
        except (ImportError, ModuleNotFoundError):
            log.warning("Failed to import translator module: %s", m)


@main.command(help="Translate metadata in supplied files and report.")
@click.argument("files", nargs=-1)
@click.option(
    "-q",
    "--quiet/--no-quiet",
    default=False,
    help="Do not report the translation content from each header. Only report failures.",
)
@hdrnum_option
@click.option(
    "-m",
    "--mode",
    default="auto",
    type=click.Choice(["auto", "verbose", "table"], case_sensitive=False),
    help="Output mode. 'verbose' prints all available information for each file found."
    " 'table' uses tabular output for a cutdown set of metadata."
    " 'auto' uses 'verbose' if one file found and 'table' if more than one is found.",
)
@regex_option
@click.pass_context
def translate(
    ctx: click.Context, files: Sequence[str], quiet: bool, hdrnum: int, mode: str, regex: str
) -> None:
    """Translate a header."""
    # For quiet mode we want to translate everything but report nothing.
    if quiet:
        mode = "none"

    okay, failed = translate_or_dump_headers(files, regex, hdrnum, ctx.obj["TRACEBACK"], output_mode=mode)

    if failed:
        click.echo("Files with failed translations:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)


@main.command(help="Dump data header to standard out in YAML format.")
@click.argument("files", nargs=-1)
@hdrnum_option
@click.option(
    "-m",
    "--mode",
    default="yaml",
    type=click.Choice(["yaml", "fixed", "yamlnative", "fixexnative"], case_sensitive=False),
    help="Output mode. 'yaml' dumps the header in YAML format (this is the default)."
    " 'fixed' dumps the header in YAML format after applying header corrections."
    " 'yamlnative' is as for 'yaml' but dumps the native (astropy vs PropertyList) native form."
    " 'yamlfixed' is as for 'fixed' but dumps the native (astropy vs PropertyList) native form.",
)
@regex_option
@click.pass_context
def dump(ctx: click.Context, files: Sequence[str], hdrnum: int, mode: str, regex: str) -> None:
    """Dump a header."""
    okay, failed = translate_or_dump_headers(files, regex, hdrnum, ctx.obj["TRACEBACK"], output_mode=mode)

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)


@main.command(help="Write JSON sidecar files alongside each data file.")
@click.argument("files", nargs=-1)
@hdrnum_option
@regex_option
@content_option
@click.pass_context
def write_sidecar(ctx: click.Context, files: Sequence[str], hdrnum: int, regex: str, content: str) -> None:
    """Write a sidecar file with header information."""
    okay, failed = write_sidecar_files(files, regex, hdrnum, content, ctx.obj["TRACEBACK"])

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay and not failed:
        # No files found at all.
        click.echo("Found no files matching regex.")
        raise click.exceptions.Exit(1)

    if not okay:
        # Good status if anything was returned in okay
        click.echo(f"No files processed successfully. Found {len(failed)}.", err=True)
        raise click.exceptions.Exit(1)


@main.command(help="Write JSON index file for entire directory.")
@click.argument("files", nargs=-1)
@hdrnum_option
@regex_option
@content_option
@click.option(
    "-o",
    "--outpath",
    type=str,
    default=None,
    help="If given, write a single index with all information in specified location."
    " Default is to write one index per directory where files are located.",
)
@click.pass_context
def write_index(
    ctx: click.Context, files: Sequence[str], hdrnum: int, regex: str, content: str, outpath: str
) -> None:
    """Write a header index file."""
    okay, failed = write_index_files(
        files, regex, hdrnum, ctx.obj["TRACEBACK"], content_mode=content, outpath=outpath
    )

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)
