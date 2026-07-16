# Glossary

This glossary provides definitions for terms used in the context of the LR Reduction software and related documentation.

## Definitions

- **run (partial)**: A run (also called a partial) is an individual measurement or data collection event within a sequence. Each run has its own unique identifier and may have specific settings or parameters. Partial data files are typically named
  `REFL_{sequence_id}_{sequence_number}_{run_number}_partial.ort`

- **sequence**: A sequence is a collection of runs that are related to a specific experiment or measurement. It represents a series of measurements taken of a particular sample at varying angles. .
  - **sequence ID**: The run number of the first run in a sequence is used as the sequence ID. This ID is used to group runs together that belong to the same sequence.
  - **sequence number**: An incremental number that identifies sequences within a set of measurements.
