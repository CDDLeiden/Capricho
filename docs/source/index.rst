CAPRICHO Documentation
======================

.. image:: ../../logo.svg
   :width: 240
   :align: center

**The ChEMBL data curator that flags issues instead of silently dropping them.**

CAPRICHO (**C**\hEMBL **A**\ggregation **P**\ackage with **R**\obust **I**\nspection and **C**\uration **H**\andling **O**\ptions) is a Python package that streamlines fetching, curating, and aggregating ChEMBL data into a machine learning-ready format for drug discovery in a flexible and reproducible manner.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   cli-reference
   guides/index
   api/index
   concepts

Goals
-----

The development of CAPRICHO is guided by two core principles:

* **Transparency Above All**: Data curation should never be a black box. Removed data points should be saved to be scrutinized by the user and the original data should be always preserved to ensure data integrity.
* **Flexibility by Design**: Every modeling project is unique. The tool must support flexible data collection and aggregation, allowing the incorporation of any ChEMBL metadata column to be incorporated into same-compound bioactivity values.

Features
--------

* Data retrieval by any ChEMBL identifier (molecule IDs, target IDs, assay IDs, or document IDs)
* Automated pChEMBL (pXC50) value calculation for bioactivities if not provided through ChEMBL
* ADMET data support with unit conversion and non-pChEMBL aggregation
* Customizable filtering options
* Configurable data aggregation options
* Save a fetching and processing recipe for reproducibility
* Command-line interface for easy use

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`