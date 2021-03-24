#!/usr/bin/env python
import os

import pytest
# our module(s)
from brukertiff_processor import BRUKERTIFFProcessor, BRUKERTIFFFile
from javabridge.jutil import JavaException
from moto import mock_dynamodb2
from moto import mock_s3
from moto import mock_ssm

from base_processor.tests import init_ssm
from base_processor.tests import setup_processor

# Test file-based exceptions
test_file_exceptions_data = [
    (
        fn,
        expected,
        {
            'file': os.path.join('/test-resources', fn)
        }
    )
    for fn, expected in [
        ('not_a_real_bruker_tiff.brukertiff.gz', JavaException),
        # ('no_xml_file', ValueError),
    ]
]


@mock_s3
@mock_dynamodb2
@pytest.mark.parametrize(
    "filename,expected,inputs",
    test_file_exceptions_data,
    ids=[
        'Java exception',
        # 'Value Error'
    ]
)
def test_file_exception_handling(filename, expected, inputs):
    with pytest.raises(expected):
        print "~" * 60
        print " Using test file %s to check exception handling when dealing with invalid NIfTI-1 files " % filename
        print "~" * 60

        mock_dynamodb2().start()
        mock_ssm().start()
        mock_s3().start()

        init_ssm()

        # init task
        tsk = BRUKERTIFFProcessor(inputs=inputs)

        setup_processor(tsk)
        print "INPUTS = ", inputs

        # Run test
        tsk.series.append(BRUKERTIFFFile())
        tsk.series[0].load_and_save_assets([], filename, prairie_version="5.2", output_filename="blah.ome.tiff")
