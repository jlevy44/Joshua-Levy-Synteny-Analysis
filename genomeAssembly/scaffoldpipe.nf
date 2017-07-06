#!/usr/bin/env nextflow

writeSh = 0

process writeShFiles {

    """
    #!/bin/bash
    python $PWD/writeShFiles.py
    ""
}
