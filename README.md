# A Geometric Clustering Tool (AGCT)
This is the Code and Dataset for AGCT and paper  "A Geometric Clustering Tool (AGCT) to robustly unravel the inner cluster structures of time-series gene expressions" by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano

## Content
Two different versions of the same program are provided. Below you can find the major differences between version 1 and version 2. You should use **version 2**.

### Version 2
 - Works with Java 8.
 - Supports highlighting sentinels (probes that are known to be significantly correlated).
 - Faster manifold and clustering computations.
 - Easier to use and better initial values.
 - More stable than version 1.

### Version 1
 - Requires Java 6 or older. Version 1 will **not** run with Java 7 or newer.
 - Allows computing the local manifold dimension estimation. This feature was removed in version 2.

The file formats are incompatible between different version. You should only use version 1 if you need to run datasets for the older version or want to see the local manifold dimension estimation.

## License

See the [LICENSE](./LICENSE) file for details

---

Thank you.