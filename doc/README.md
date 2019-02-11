# How to build the documentation

This documentation is built with [Documenteer](https://documenteer.lsst.io), which wraps [Sphinx](http://www.sphinx-doc.org/en/stable/). To install the documentation build dependencies:

```sh
pip install -r requirements.txt
```

[Build the documentation with Documenteer](https://documenteer.lsst.io/pipelines/package-docs-cli.html):

```sh
package-docs build
```

If necessary, you can also clean the built documentation and intermediate products:

```sh
package-docs clean
```

Travis CI builds the documentation and deploys it to an appropriate edition on LSST the Docs. See https://astro-metadata-translator.lsst.io/v for available editions.
