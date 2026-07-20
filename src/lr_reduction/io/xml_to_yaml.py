"""Convert a legacy RefRed XML reduction template into the new-workflow YAML configuration format (§3.3.10).

This is the sanctioned way to bring a legacy XML template into the new
workflow: the new workflow's own loader (`lr_reduction.io.config_loader.ConfigLoader`) never reads
XML directly for that purpose (`io.xml.read_config` is a separate, deprecated backward-compatibility
path). Run this utility first, then pass the resulting YAML file to the new workflow.

Available both as a library call and as a command-line entry point — the installed
`lr-xml-to-yaml` console script, or `python -m lr_reduction.io.xml_to_yaml` from a repository
checkout without installing the package first:

    >>> from lr_reduction.io.xml_to_yaml import convert
    >>> convert("template.xml", "config.yaml")
"""

import argparse

import yaml

from lr_reduction.io.xml import read_config


def convert(input_path: str, output_path: str) -> None:
    """Read a legacy XML template from *input_path* and write it as YAML to *output_path*."""
    config = read_config(input_path)
    with open(output_path, "w") as f:
        yaml.safe_dump(config.model_dump(exclude_none=True), f, sort_keys=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to the legacy XML reduction template")
    parser.add_argument("output", help="Path to write the converted YAML configuration to")
    args = parser.parse_args()
    convert(args.input, args.output)


if __name__ == "__main__":
    main()
