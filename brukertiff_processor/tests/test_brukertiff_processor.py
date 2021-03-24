#!/usr/bin/env python
import os
import glob

import pytest
# our module(s)
from brukertiff_processor import BRUKERTIFFProcessor
from moto import mock_s3
from moto import mock_ssm

from base_processor.tests import init_ssm
from base_processor.tests import setup_processor

test_processor_data = [
    'TSeries-06272017-1624-103.brukertiff.gz'
]


@pytest.mark.parametrize("dirname", test_processor_data)
def test_brukertiff_processor(dirname):
    mock_ssm().start()
    mock_s3().start()

    init_ssm()

    # init task
    inputs = {'file': os.path.join('/test-resources', dirname)}
    task = BRUKERTIFFProcessor(inputs=inputs)

    setup_processor(task)

    # run
    task.run()

    # Confirm output files
    assert len(task.series) == 3
    assert len(glob.glob('*.ome.tiff')) == 3
    assert os.path.isfile('metadata.json')

    # Confirm view asset json
    assert len(glob.glob('*-zoomed')) == 3

    # Clean up
    os.system('rm *zoomed/*.png')
    os.system('rm *zoomed/*.dzi')

    mock_s3().stop()
    mock_ssm().stop()
