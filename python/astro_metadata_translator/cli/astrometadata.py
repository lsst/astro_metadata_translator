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

__all__ = ("main",)

import importlib
import logging
import click

from ..bin.translateheader import process_files as translate_header
from ..bin.writesidecar import write_sidecar_files
from ..bin.writeindex import write_index_files

# Default regex for finding data files
re_default = r"\.fit[s]?\b"


@click.group(name="astrometadata", context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--log-level",
              type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False),
              default="INFO",
              help="Python logging level to use.")
@click.option("--traceback/--no-traceback", default=False,
              help="Give detailed trace back when any errors encountered.")
@click.option("-p", "--packages", multiple=True,
              help="Python packages to import to register additional translators")
@click.pass_context
def main(ctx, log_level, traceback, packages):
    ctx.ensure_object(dict)

    logging.basicConfig(level=log_level)

    # Traceback needs to be known to subcommands
    ctx.obj["TRACEBACK"] = traceback

    # Process import requests
    for m in packages:
        importlib.import_module(m)


@main.command(help="Translate metadata in supplied files and report.")
@click.argument("files", nargs=-1)
@click.option("-q", "--quiet/--no-quiet",
              default=False,
              help="Do not report the translation content from each header. Only report failures.")
@click.option("-n", "--hdrnum", default=1,
              help="HDU number to read. If the HDU can not be found, a warning is issued but translation"
              " is attempted using the primary header. The primary header is always read and merged with"
              " this header.")
@click.option("-m", "--mode",
              default="auto",
              type=click.Choice(["auto", "verbose", "table"], case_sensitive=False),
              help="Output mode. 'verbose' prints all available information for each file found."
              " 'table' uses tabular output for a cutdown set of metadata."
              " 'auto' uses 'verbose' if one file found and 'table' if more than one is found.")
@click.option("-r", "--regex",
              default=re_default,
              help="When looking in a directory, regular expression to use to determine whether"
              f" a file should be examined. Default: '{re_default}'")
@click.pass_context
def translate(ctx, files, quiet, hdrnum, mode, regex):

    # For quiet mode we want to translate everything but report nothing.
    if quiet:
        mode = "none"

    okay, failed = translate_header(files, regex, hdrnum, ctx.obj["TRACEBACK"], output_mode=mode)

    if failed:
        click.echo("Files with failed translations:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)


@main.command(help="Dump data header to standard out in YAML format.")
@click.argument("files", nargs=-1)
@click.option("-n", "--hdrnum", default=1,
              help="HDU number to read. If the HDU can not be found, a warning is issued but translation"
              " is attempted using the primary header. The primary header is always read and merged with"
              " this header.")
@click.option("-m", "--mode",
              default="yaml",
              type=click.Choice(["yaml", "fixed", "yamlnative", "fixexnative"], case_sensitive=False),
              help="Output mode. 'yaml' dumps the header in YAML format (this is the default)."
              " 'fixed' dumps the header in YAML format after applying header corrections."
              " 'yamlnative' is as for 'yaml' but dumps the native (astropy vs PropertyList) native form."
              " 'yamlfixed' is as for 'fixed' but dumps the native (astropy vs PropertyList) native form.")
@click.option("-r", "--regex",
              default=re_default,
              help="When looking in a directory, regular expression to use to determine whether"
              f" a file should be examined. Default: '{re_default}'")
@click.pass_context
def dump(ctx, files, hdrnum, mode, regex):

    okay, failed = translate_header(files, regex, hdrnum, ctx.obj["TRACEBACK"], output_mode=mode)

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)


@main.command(help="Write JSON sidecar files with ObservationInfo summary information.")
@click.argument("files", nargs=-1)
@click.option("-n", "--hdrnum", default=1,
              help="HDU number to read. If the HDU can not be found, a warning is issued but translation"
              " is attempted using the primary header. The primary header is always read and merged with"
              " this header.")
@click.option("-r", "--regex",
              default=re_default,
              help="When looking in a directory, regular expression to use to determine whether"
              f" a file should be examined. Default: '{re_default}'")
@click.pass_context
def write_sidecar(ctx, files, hdrnum, regex):
    okay, failed = write_sidecar_files(files, regex, hdrnum, ctx.obj["TRACEBACK"])

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)


@main.command(help="Write JSON index file for entire directory.")
@click.argument("files", nargs=-1)
@click.option("-n", "--hdrnum", default=1,
              help="HDU number to read. If the HDU can not be found, a warning is issued but translation"
              " is attempted using the primary header. The primary header is always read and merged with"
              " this header.")
@click.option("-r", "--regex",
              default=re_default,
              help="When looking in a directory, regular expression to use to determine whether"
              f" a file should be examined. Default: '{re_default}'")
@click.pass_context
def write_index(ctx, files, hdrnum, regex):
    okay, failed = write_index_files(files, regex, hdrnum, ctx.obj["TRACEBACK"])

    if failed:
        click.echo("Files with failed header extraction:", err=True)
        for f in failed:
            click.echo(f"\t{f}", err=True)

    if not okay:
        # Good status if anything was returned in okay
        raise click.exceptions.Exit(1)
