import fcntl
import itertools as it
import json
import os
import re
import sys
import traceback
from distutils.version import LooseVersion
from xml.etree import ElementTree

import PIL.Image
import bioformats as bf
import bioformats.omexml as ome
import numpy as np
from base_processor.imaging import utils
from javabridge.jutil import JavaException

from base_image_microscopy_processor import BaseMicroscopyImageProcessor
from utils import get_prairieview_version, parse_cfg_file, get_element_size_um, timestamp

reload(sys)
sys.setdefaultencoding('utf8')


class BRUKERTIFFFile(object):
    def __init__(self, *args, **kwargs):
        """ Constructor for Bruker PrairieView TIFF image sets """
        self.img_dimensions = {}
        self.num_dimensions = -1
        self.file_path = None
        self.view_path = None
        self.ometiff_path = None

        self.img_data = kwargs.get('img_data', None)
        self.img_data_dtype = kwargs.get('img_data_dtype', None)
        self.metadata = kwargs.get('metadata', {})

        self.img_rdr = kwargs.get('img_rdr', None)
        self.DimensionOrder = kwargs.get('DimensionOrder', None)
        self.ImageCount = kwargs.get('ImageCount', -1)
        self.SizeX = kwargs.get('SizeX', -1)
        self.SizeY = kwargs.get('SizeY', -1)
        self.SizeZ = kwargs.get('SizeZ', -1)
        self.SizeC = kwargs.get('SizeC', -1)
        self.SizeT = kwargs.get('SizeT', -1)
        self.PixelType = kwargs.get('PixelType', None)
        self.RGBChannelCount = kwargs.get('RGBChannelCount', -1)

        self.isRGB = kwargs.get('isRGB', False)
        self.RGBDimension = kwargs.get('RGBDimension', -1)
        self.hasTimeDimension = kwargs.get('hasTimeDimension', False)
        self.TimeDimension = kwargs.get('TimeDimension', -1)

        self.view_format = 'dzi'
        self.optimize = kwargs.get('optimize', False)
        self.tile_size = kwargs.get('tile_size', 128)
        self.tile_overlap = kwargs.get('tile_overlap', 0)
        self.tile_format = kwargs.get('tile_format', "png")
        self.image_quality = kwargs.get('image_quality', 1.0)
        self.resize_filter = kwargs.get('resize_filter', "bicubic")

    def _load_and_save_assets(self, i_xyzct, n_xyzct, tiffs):
        """
        Loads image, populated image data matrix and sets image properties using bioformats ImageReader.

        Args:
            tiffs: Multi-dimensional list of tiff files in directory.
            xml_file_path: Absolute path to config filename.
            prairie_version: PrairieView prairie_version number.

        Returns:

        Raises:
            JavaException: Raises a Java error exception when bioformats encounters an error reading the config XML
            file or OME-TIFF files in the PrairiveView Image set directory.
            TypeError: Raises an XML decoding error (for e.g. encoding issues)
        """

        # Parse parallelization arguments
        i_x, i_y, i_z, i_c, i_t = i_xyzct
        n_x, n_y, n_z, n_c, n_t = n_xyzct

        if i_x == -1 or n_x == -1:
            i_x = 0
            n_x = 1
        if i_y == -1 or n_y == -1:
            i_y = 0
            n_y = 1
        if i_z == -1 or n_z == -1:
            i_z = 0
            n_z = 1
        if i_c == -1 or n_c == -1:
            i_c = 0
            n_c = 1
        if i_t == -1 or n_t == -1:
            i_t = 0
            n_t = 1

        # Make view directory
        if not os.path.exists('%s-zoomed' % os.path.basename(self.file_path)):
            os.makedirs('%s-zoomed' % os.path.basename(self.file_path))
        asset_format = self.view_format

        for z in range(int(self.SizeZ / float(n_z) * i_z), int(self.SizeZ / float(n_z) * (i_z + 1))):
            for t, tiff in enumerate(tiffs):
                for c in range(len(tiff[0])):
                    if t not in range(int(self.SizeT / float(n_t) * i_t),
                                      int(self.SizeT / float(n_t) * (i_t + 1))):
                        continue
                    if c not in range(int(self.SizeC / float(n_c) * i_c),
                                      int(self.SizeC / float(n_c) * (i_c + 1))):
                        continue

                    # Get image data matrix
                    # Y and X flipped because image
                    self.img_data = np.zeros((self.SizeY, self.SizeX))
                    with bf.ImageReader(path=tiff[z][c]) as rdr:
                        self.img_data = rdr.read(z=z)

                    # Convert type of image in order to save in deep-zoom
                    self.img_data = self._convert_image_data_type(img_data=self.img_data)
                    image = PIL.Image.fromarray(self.img_data)
                    self._save_view(image, z=z, c=c, t=t, is_rgb=True, asset_format=asset_format)

        self.view_path = os.path.join(os.getcwd(), '%s-zoomed' % os.path.basename(self.file_path))
        return

    def _convert_image_data_type(self, img_data=None, format=np.uint8):
        """Convert image data matrix datatype to another format"""

        # Check correct data type format to convert
        if format in np.typeDict.values():
            pass
        elif type(format) == str:
            try:
                format = np.dtype(format)
            except Exception:
                raise TypeError("Image data type to convert to should be of type numpy.dtype")
        else:
            raise TypeError("Image data type to convert to should be of type numpy.dtype")

        if img_data is None:
            # Convert only if image data type is not in desired data type format
            if self.img_data_dtype != format:
                self.img_data = utils.convert_image_data_type(self.img_data, format)
                self.img_data_dtype = format
            return
        else:
            img_data = utils.convert_image_data_type(img_data, format)
            return img_data

    def _save_view(self, image, z=None, c=None, t=None, is_rgb=False, asset_format='dzi'):
        """"""

        if is_rgb:
            # Save asset in appropriate format
            filename = os.path.join(
                '%s-zoomed' % os.path.basename(self.file_path),
                'dim_Z_slice_{slice_z}_dim_T_slice_{slice_t}.{fmt}'.format(
                    slice_z=z, slice_t=t, fmt=asset_format)
            )

            utils.save_asset(
                image,
                asset_format,
                filename,
                optimize=self.optimize, tile_size=self.tile_size,
                tile_overlap=self.tile_overlap, tile_format=self.tile_format,
                image_quality=self.image_quality, resize_filter=self.resize_filter
            )
            # Save thumbnail
            if asset_format == 'dzi':
                timage = image.copy()
                timage.thumbnail((200, 200), PIL.Image.ANTIALIAS)
                timage.save(
                    os.path.join(
                        '%s-zoomed' % os.path.basename(self.file_path),
                        'dim_Z_slice_{slice_z}_dim_T_slice_{slice_t}_thumbnail.png'.format(
                            slice_z=z,
                            slice_t=t
                        )
                    )
                )
                # Create large thumbnail
                timage = image.copy()
                timage.thumbnail((1000, 1000), PIL.Image.ANTIALIAS)
                timage.save(
                    os.path.join(
                        '%s-zoomed' % os.path.basename(self.file_path),
                        'dim_Z_slice_{slice_z}_dim_T_slice_{slice_t}_large_thumbnail.png'.format(
                            slice_z=z,
                            slice_t=t
                        )
                    )
                )
        else:
            # Save asset in appropriate format
            filename = os.path.join(
                '%s-zoomed' % os.path.basename(self.file_path),
                'dim_Z_slice_{slice_z}_dim_C_slice_{slice_c}_dim_T_slice_{slice_t}.{fmt}'.format(
                    slice_z=z, slice_c=c, slice_t=t, fmt=asset_format)
            )

            utils.save_asset(
                image,
                asset_format,
                filename,
                optimize=self.optimize, tile_size=self.tile_size,
                tile_overlap=self.tile_overlap, tile_format=self.tile_format,
                image_quality=self.image_quality, resize_filter=self.resize_filter
            )
            # Save thumbnail
            if asset_format == 'dzi':
                timage = image.copy()
                timage.thumbnail((200, 200), PIL.Image.ANTIALIAS)
                timage.save(
                    os.path.join(
                        '%s-zoomed' % os.path.basename(self.file_path),
                        'dim_Z_slice_{slice_z}_dim_C_slice_{slice_c}_'
                        'dim_T_slice_{slice_t}_thumbnail.png'.format(
                            slice_z=z,
                            slice_c=c,
                            slice_t=t
                        )
                    )
                )
                # Create large thumbnail
                timage = image.copy()
                timage.thumbnail((1000, 1000), PIL.Image.ANTIALIAS)
                timage.save(
                    os.path.join(
                        '%s-zoomed' % os.path.basename(self.file_path),
                        'dim_Z_slice_{slice_z}_dim_C_slice_{slice_c}_'
                        'dim_T_slice_{slice_t}_large_thumbnail.png'.format(
                            slice_z=z,
                            slice_c=c,
                            slice_t=t
                        )
                    )
                )
        return

    def set_img_properties(self, tiffs, xml_file_path, prairie_version='5.2', asset_format='dzi'):
        """Create and assign properties for the image"""

        self.view_format = asset_format

        # Load image using BioFormats
        try:
            ImageReader = bf.formatreader.make_image_reader_class()
            self.img_rdr = ImageReader()
            self.img_rdr.setId(xml_file_path)
        except JavaException:
            raise

        # Get image data details
        self.DimensionOrder = list(self.img_rdr.getDimensionOrder())
        self.ImageCount = self.img_rdr.getImageCount()
        self.SizeX = self.img_rdr.getSizeX()
        self.SizeY = self.img_rdr.getSizeY()
        self.SizeC = self.img_rdr.getSizeC()
        self.SizeT = len(tiffs)
        self.SizeZ = len(tiffs[0])

        self.PixelType = self.img_rdr.getPixelType()
        self.RGBChannelCount = self.img_rdr.getRGBChannelCount()
        self.isRGB = self.img_rdr.isRGB()

        self.RGBDimension = 3
        self.hasTimeDimension = True
        self.TimeDimension = 4

        # Set data type #TODO: For some reason, rdr.read returns np.float64 instead of specified PixelType
        self.PixelType = 6
        if self.PixelType == 0:  # int8
            # self.img_data = self.img_data.astype(np.int8)
            self.img_data_dtype = np.int8
        elif self.PixelType == 1:  # uint8
            # self.img_data = self.img_data.astype(np.uint8)
            self.img_data_dtype = np.uint8
        elif self.PixelType == 2:  # int16
            # self.img_data = self.img_data.astype(np.int16)
            self.img_data_dtype = np.int16
        elif self.PixelType == 3:  # uint16
            # self.img_data = self.img_data.astype(np.uint16)
            self.img_data_dtype = np.uint16
        elif self.PixelType == 4:  # int32
            # self.img_data = self.img_data.astype(np.int32)
            self.img_data_dtype = np.int32
        elif self.PixelType == 5:  # uint32
            # self.img_data = self.img_data.astype(np.uint32)
            self.img_data_dtype = np.uint32
        elif self.PixelType == 6:  # float
            # self.img_data = self.img_data.astype(np.float)
            self.img_data_dtype = np.float
        elif self.PixelType == 7:  # double
            # self.img_data = self.img_data.astype(np.double)
            self.img_data_dtype = np.double

        # Set number of dimensions of image matrix
        self.num_dimensions = 5
        self.img_data_shape = [self.SizeX, self.SizeY, self.SizeZ, self.SizeC, self.SizeT]

        dim_assignment = list('XYZCT')  # Force assignment

        self.img_dimensions['filename'] = os.path.basename(self.file_path)
        self.img_dimensions['num_dimensions'] = self.num_dimensions
        self.img_dimensions['isColorImage'] = False
        self.img_dimensions['dimensions'] = {}

        for dim in range(self.num_dimensions):
            self.img_dimensions['dimensions'][dim] = {}
            self.img_dimensions['dimensions'][dim]["assignment"] = dim_assignment[dim]
            self.img_dimensions['dimensions'][dim]["length"] = self.img_data_shape[dim]
            self.img_dimensions['dimensions'][dim]["resolution"] = -1
            self.img_dimensions['dimensions'][dim]["units"] = "um"
            if dim_assignment[dim] == 'C' and self.isRGB:
                self.RGBDimension = dim
            if dim_assignment[dim] == 'T':
                self.hasTimeDimension = True
                self.TimeDimension = dim
        self.img_dimensions['isColorImage'] = self.isRGB

        # Get metadata
        metadata = bf.get_omexml_metadata(xml_file_path).encode('utf-8')

        # Get resolution
        # Get the size in um of the x and y dimensions
        element_size_um = get_element_size_um(
            xml_file_path, prairie_version)
        self.img_dimensions['dimensions'][0]['resolution'] = element_size_um[1]  # X dim
        self.img_dimensions['dimensions'][0]['resolution_units'] = 'um'
        self.img_dimensions['dimensions'][1]['resolution'] = element_size_um[2]  # X dim
        self.img_dimensions['dimensions'][1]['resolution_units'] = 'um'


        resx_node = 'PhysicalSizeX="'
        resy_node = 'PhysicalSizeY="'

        start = metadata.index(resx_node) + len(resx_node)
        end = metadata.index('"', start)
        assert self.img_dimensions['dimensions'][0]["resolution"] == float(metadata[start:end])

        start = metadata.index(resy_node) + len(resy_node)
        end = metadata.index('"', start)
        assert self.img_dimensions['dimensions'][1]["resolution"] == float(metadata[start:end])

        try:
            self.metadata = bf.omexml.OMEXML(metadata.decode('utf-8'))
        except Exception:
            # Raise an exception only if we file ends in .ome.tif/.ome.tiff since we know it should be OME compatible
            if self.file_path.endswith('.czi'):
                raise TypeError("XML metadata decoding failed.")
        # self.metadata = javabridge.jdictionary_to_string_dictionary(self.img_rdr.getMetadata())
        return

    def load_and_save_assets(self, tiffs, xml_path, prairie_version, output_filename,
                             i_xyzct=(-1, -1, -1, -1, -1), n_xyzct=(-1, -1, -1, -1, -1), asset_format='dzi'):
        """Load image and generate view assets"""
        # Set file path
        self.file_path = output_filename

        # Set image properties
        self.set_img_properties(tiffs=tiffs, xml_file_path=xml_path, prairie_version=prairie_version, asset_format=asset_format)

        # Load image
        self._load_and_save_assets(i_xyzct=i_xyzct, n_xyzct=n_xyzct, tiffs=tiffs)

    def save_ometiff(self, output_filename, tiffs):
        """Save OME-TIFF images derived from CZI file"""

        if self.img_data_dtype == np.int8:
            data_type = ome.PT_INT8
        elif self.img_data_dtype == np.uint8:
            data_type = ome.PT_UINT8
        elif self.img_data_dtype == np.int16:
            data_type = ome.PT_INT16
        elif self.img_data_dtype == np.uint16:
            data_type = ome.PT_UINT16
        elif self.img_data_dtype == np.int32:
            data_type = ome.PT_INT32
        elif self.img_data_dtype == np.uint32:
            data_type = ome.PT_UINT32
        elif self.img_data_dtype == np.float:
            data_type = ome.PT_FLOAT
        elif self.img_data_dtype == np.double:
            data_type = ome.PT_DOUBLE
        else:
            raise NotImplementedError

        for z in range(self.SizeZ):
            for t, tiff in enumerate(tiffs):
                for c in range(len(tiff[0])):

                    # Get image data matrix
                    # Y and X flipped because image
                    self.img_data = np.zeros((self.SizeY, self.SizeX))
                    with bf.ImageReader(path=tiff[z][c]) as rdr:
                        self.img_data = rdr.read(z=z)

                    # Convert type of image in order to save in deep-zoom
                    self.img_data = self._convert_image_data_type(img_data=self.img_data)

                    # Write OMETIFF
                    bf.formatwriter.write_image(output_filename, self.img_data, data_type,
                                                z=z, c=c, t=t,
                                                size_c=self.SizeC, size_z=self.SizeZ, size_t=self.SizeT)
        self.ometiff_path = output_filename
        return

    @property
    def view_size(self):
        """ Returns the full size of view assets directory """
        return os.path.getsize(self.view_path)

    def get_view_asset_dict(self, storage_bucket, upload_key):
        """ Returns the view mutation as JSON object """
        upload_key = upload_key.rstrip('/')
        json_dict = {
            "bucket": storage_bucket,
            "key": upload_key,
            "type": "View",
            "size": self.view_size,
            "fileType": "Image"
        }
        return json_dict

    def get_dim_assignment(self):
        """ Return dimension assignment for each dimension as a list of str. """

        return list(self.DimensionOrder)




class BRUKERTIFFProcessor(BaseMicroscopyImageProcessor):
    required_inputs = ['file']

    def __init__(self, *args, **kwargs):
        """Constructor for Bruker PrairieView TIFF Processor"""

        super(BRUKERTIFFProcessor, self).__init__(*args, **kwargs)
        # Get input "file", which should be a imaging study directory
        self.file = self.inputs.get('file')
        self.unzipped_file_dir = None

        # Initialize all the series of images that the processor will decode from study directory
        self.series = []

        # Initialize image asset properties
        self.upload_key = None
        try:
            self.optimize = utils.str2bool(self.inputs.get('optimize_view'))
        except AttributeError:
            self.optimize = False

        try:
            self.tile_size = int(self.inputs.get('tile_size'))
        except (ValueError, KeyError, TypeError) as e:
            self.tile_size = 128

        try:
            self.tile_overlap = int(self.inputs.get('tile_overlap'))
        except (ValueError, KeyError, TypeError) as e:
            self.tile_overlap = 0

        try:
            self.tile_format = self.inputs.get('tile_format')
            if self.tile_format is None:
                self.tile_format = "png"
        except KeyError:
            self.tile_format = "png"

        try:
            self.image_quality = float(self.inputs.get('image_quality'))
        except (ValueError, KeyError, TypeError) as e:
            self.image_quality = 1.0

        try:
            self.resize_filter = self.inputs.get('resize_filter')
        except KeyError:
            self.resize_filter = "bicubic"

    def parse_prairieview_xml(self, xml_path):
        """
        Parse Bruker PrairieView XML header file to get PrairieView version, config file name and XML elements

        Args:
            xml_path: File path to XML header file

        Returns:
            prairie_version: version (str)
            cfg_filename: Config filename (str)
            protocol_elements: Different protocol elements parsed from XML header file (list)
            cycles: Different image acquisition cycles parsed from XML header file (list)

        Raises:
            IOError: raised if get_prairieview_version fails on decoding XML file
            ValueError: if Prairie version is too old (not supported)
        """
        xml_file_name = os.path.basename(xml_path)
        image_directory = os.path.dirname(xml_path)
        try:
            prairie_version = get_prairieview_version(xml_path)
        except (IOError, ElementTree.ParseError):
            self.LOGGER.error("{} XML Parse error: {}".format(timestamp(), image_directory))
            self.LOGGER.error(traceback.format_exc())
            raise IOError

        # If data was recorded pre-5.0, just skip the folder
        if prairie_version < LooseVersion('5.0'):
            failed_message = "{} Prairie version too old ({}): {}".format(
                timestamp(), prairie_version, image_directory)
            self.LOGGER.error(failed_message)
            raise ValueError(failed_message)

        if prairie_version > LooseVersion('5.2'):
            cfg_filename = xml_file_name.replace('.xml', '.env')
        else:
            cfg_filename = xml_file_name.replace('.xml', 'Config.cfg')
        self.LOGGER.info('Since Prairie version is %s, using config file %s' % (prairie_version, cfg_filename))

        # Obtain different elements and number of iterations
        protocol_elements, n_cycles = parse_cfg_file(
            os.path.join(image_directory, cfg_filename))

        # Older Prairie versions don't actually store iterations/cycles
        # anywhere...hack to not need to know beforehand
        if n_cycles is None:
            cycles = it.count()
        else:
            cycles = range(n_cycles)

        return prairie_version, cfg_filename, protocol_elements, cycles

    def collect_tiffs_to_save(self, cycles, protocol_elements, sequences, prairie_version, xml_path):
        """
        Collect all the tiffs inside of a image set directory that will be decoded and stacked

        Args:
            cycles: Different image acquisition cycles parsed from XML header file (list)
            protocol_elements: Different protocol elements parsed from XML header file (list)
            sequences:
            prairie_version: PrairieView version based on XML header file
            xml_path: Path to XML header file

        Returns:
            tiffs_to_save: Dictionary with
                key= output OME-TIFF filename (str) and value= list of tiffs and cycles (tuple)
            failed_message: Any error message raised during collection of all TIFFs to stack into output file

        Raises:

        """
        image_directory = os.path.dirname(xml_path)
        basename = os.path.basename(xml_path).replace('.xml', '')
        failed_message = False
        tiffs_to_save = {}
        for cycle in cycles:
            for idx, (protocol, reps) in enumerate(protocol_elements):
                output_filename = '{}_Cycle{:05d}_Element{:05d}.ome.tiff'.format(basename, cycle + 1, idx + 1)

                self.LOGGER.info('Writing to output filename {}'.format(output_filename))

                # Determine tiff files that are under protocol 'PVTSeriesElementSequence'
                if protocol == 'PVTSeriesElementSequence':
                    # This will both check for mis-matched number of
                    # sequences and fix for old Prairie not knowing cycles
                    try:
                        sequence = sequences.next()
                    except StopIteration:
                        if prairie_version < LooseVersion('5.2'):
                            failed_message = True
                            self.LOGGER.warn(
                                'Stop Iteration exception raised in Prairie version %s. '
                                'Stopping sequence iteration.' % prairie_version)
                            break
                        else:
                            err_msg = '{} Sequence length mis-match, '.format(
                                timestamp()) + \
                                      '{} expected, {} actual: {}'.format(
                                          len(protocol_elements), idx + 1, image_directory)
                            self.LOGGER.error(err_msg)
                            failed_message = err_msg
                            break
                    frames = sequence.findall('Frame')
                    channels = [ff.get('channelName')
                                for ff in frames[0].findall('File')]
                    channels.sort()
                    if len(frames) != reps:
                        err_msg = '{} Frame/rep '.format(timestamp()) \
                                  + 'mismatch, {} frames, {} reps: {}'.format(
                            len(frames), reps, image_directory)
                        self.LOGGER.error(err_msg)
                        # TODO: Deal with mismatch
                        self.LOGGER.info('Ignoring mismatch for now.')
                    tiff_files = []
                    # Each frame is a time step
                    for frame in frames:
                        files = [os.path.join(image_directory, ff.get('filename'))
                                 for ff in frame.findall('File')]
                        files.sort()
                        tiff_files.append([files])

                # Determine tiff files that are under protocol 'PVTSeriesElementZSeries'
                elif protocol == 'PVTSeriesElementZSeries':
                    tiff_files = []
                    # Each sequence/rep is a time step
                    channels = None
                    for rep in range(reps):
                        try:
                            sequence = sequences.next()
                        except StopIteration:
                            # If we run out of sequences on the first rep,
                            # there's just no more cycles left, which is fine.
                            # If happens in the middle of reps we are actually
                            # missing data.
                            if prairie_version < LooseVersion('5.2') and rep == 0:
                                failed_message = True
                                self.LOGGER.warn(
                                    'Stop Iteration exception raised in Prairie version %s. '
                                    'Stopping sequence iteration.' % prairie_version)
                                break
                            else:
                                err_msg = \
                                    '{} Sequence length mis-match, '.format(
                                        timestamp()) \
                                    + '{} expected, {} actual: {}'.format(
                                        reps, rep + 1, image_directory)
                                self.LOGGER.error(err_msg)
                                failed_message = err_msg
                                break
                        frames = sequence.findall('Frame')
                        if channels is None:
                            channels = [ff.get('channelName')
                                        for ff in frames[0].findall('File')]
                            channels.sort()
                        tiff_files.append([])
                        # Each frame is a z-plane
                        for frame in frames:
                            files = [os.path.join(image_directory, ff.get('filename'))
                                     for ff in frame.findall('File')]
                            files.sort()
                            tiff_files[-1].append(files)
                    if failed_message:
                        break

                # Unrecognized protocol element
                else:
                    err_msg = '{} Unrecognized '.format(timestamp()) + \
                              'protocol element, skipping directory: {}, {}'.format(
                                  image_directory, protocol)
                    self.LOGGER.error(err_msg)
                    failed_message = err_msg
                    break

                # Add to dictionary of all tiffs to save
                tiffs_to_save[output_filename] = (tiff_files, channels)

            # Break loop if any errors: e.g. mismatch between config file and list of actual tiff files,
            # missing tiff files, et.c
            if failed_message:
                break

        return tiffs_to_save, failed_message

    def check_tiffs_to_save(self, tiffs_to_save, sequences, failed):
        """
        Check to see if all tiffs were collected correctly before saving OME-TIFF output
        files and exploding view assets.

        Args:
            tiffs_to_save: Dictionary with
                key= output OME-TIFF filename (str) and value= list of tiffs and cycles (tuple)
            sequences:
            failed: Any error message raised during collection of all TIFFs to stack into output file

        Returns:

        Raises:
            Exception: raised if any exceptions were raised during collection of tiffs or if not all files
                are present in image directory
            ValueError: raised if sequence length mismatched in sequences

        """
        # Log reported errors during iteration of cycles/protocols
        if failed:
            if len(failed) == 3:
                self.LOGGER.error(traceback.print_tb(failed[2]))
            else:
                self.LOGGER.error(failed)
            raise Exception

        # Check if there's nothing to do
        if not len(tiffs_to_save):
            raise ValueError

        # Make sure we've exactly gone through all of the sequences
        try:
            sequences.next()
        except StopIteration:
            pass
        else:
            err_msg = '{} Sequence length mis-matching'.format(timestamp()) + \
                      ', skipping directory: {}'.format(self.unzipped_file_dir)
            self.LOGGER.error(err_msg)
            raise ValueError(err_msg)

        # Check to make sure all of the files are there
        for files, _ in tiffs_to_save.itervalues():
            for frame in files:
                for z_plane in frame:
                    for f in z_plane:
                        if not os.path.exists(f):
                            err_msg = '{} Missing file'.format(timestamp()) + \
                                      ', skipping directory: {}'.format(f)
                            self.LOGGER.error(err_msg)
                            failed = err_msg
                            break
                    if failed:
                        break
                if failed:
                    break
            if failed:
                break
        if failed:
            self.LOGGER.error(failed)
            raise Exception

    def save_outputs(self, tiffs_to_save, xml_path, prairie_version, i_xyzct, n_xyzct):
        """
        Save output OME-TIFF files and associated exploded view assets

        Args:
            tiffs_to_save: Dictionary with
                key= output OME-TIFF filename (str) and value= list of tiffs and cycles (tuple)
            xml_path: File path to XML header file
            prairie_version: Version of PrairieView image set decoded from XML header

        Returns:
            failed_message: List of all exceptions raised with associated messages

        Raises:

        """
        image_directory = os.path.dirname(xml_path)

        # For all series uncovered in generation of list of tiff files:
        # 1. Add it to self.series
        # 2. Load that series image
        # 3. Convert to appropriate image data type in preparation for generation of view assets
        # 4. Save to local storage in 'view/' directory prior to publishing of outputs or upload to S3
        failed_message = False
        for series_index, (output_filename, (tiffs, channels)) in enumerate(tiffs_to_save.iteritems()):
            self.LOGGER.info("{} Creating {}".format(
                timestamp(), os.path.join(image_directory, output_filename)))
            try:
                # Initialize output file and image dimensions
                self.series.append(BRUKERTIFFFile(optimize=self.optimize, tile_size=self.tile_size,
                                                  tile_overlap=self.tile_overlap, tile_format=self.tile_format,
                                                  image_quality=self.image_quality, resize_filter=self.resize_filter))
                self.series[series_index].load_and_save_assets(
                    tiffs,
                    xml_path,
                    output_filename=output_filename,
                    prairie_version=prairie_version,
                    i_xyzct=i_xyzct,
                    n_xyzct=n_xyzct
                )

                # Create output OME-TIFF
                self.LOGGER.info('!!!! WRITING OMETIFF NOW with %s' % output_filename)

                if i_xyzct == (-1, -1, -1, -1, -1) or i_xyzct == (0, 0, 0, 0, 0):
                    self.series[series_index].save_ometiff(output_filename, tiffs)
            except Exception:
                self.LOGGER.error(("{} FAILED creating {}".format(
                    timestamp(), output_filename)))
                failed_message = sys.exc_info()
                self.LOGGER.error(failed_message)
            else:
                self.LOGGER.info("{} Successfully created {}".format(
                    timestamp(), output_filename))

                ## Hack for backwards compatibility for DeepZoom viewer
                # Create RGB image regardless of isRGB depending on number of channel count
                if i_xyzct == (-1, -1, -1, -1, -1) or i_xyzct == (0, 0, 0, 0, 0):
                    img_data = np.zeros(
                        (self.series[series_index].SizeX, self.series[series_index].SizeY, self.series[series_index].SizeC))
                    for z in range(self.series[series_index].SizeZ):
                        for t, tiff in enumerate(tiffs):
                            for c in range(len(tiff[0])):
                                # Get image data matrix
                                # Y and X flipped because image
                                self.series[series_index].img_data = np.zeros((self.series[series_index].SizeY, self.series[series_index].SizeX))
                                with bf.ImageReader(path=tiff[z][c]) as rdr:
                                    img_data_slice = rdr.read(z=z)
                                    img_data[:, :, c] = np.max((img_data[:, :, c].squeeze(), img_data_slice), axis=0)
                    mip_img_data = np.zeros((self.series[series_index].SizeX, self.series[series_index].SizeY, 3))
                    if self.series[series_index].SizeC == 1:
                        mip_img_data[:, :, :1] = img_data
                    if self.series[series_index].SizeC == 2:
                        mip_img_data[:, :, :2] = img_data
                    if self.series[series_index].SizeC == 3:
                        mip_img_data = img_data
                    if self.series[series_index].SizeC > 3:
                        mip_img_data = img_data[:, :, :3]

                    mip_img_data = utils.convert_image_data_type(mip_img_data, np.uint8)
                    image = PIL.Image.fromarray(mip_img_data, "RGB")

                    # Save asset as slide.dzi for backwards compatibility
                    filename = os.path.join(
                        '%s-zoomed' % os.path.basename(self.series[series_index].file_path),
                        'slide.dzi'
                    )
                    utils.save_asset(
                        image,
                        'dzi',
                        filename,
                        optimize=self.optimize, tile_size=self.tile_size,
                        tile_overlap=self.tile_overlap, tile_format=self.tile_format,
                        image_quality=self.image_quality, resize_filter=self.resize_filter
                    )

                    # Save thumbnail
                    timage = image.copy()
                    timage.thumbnail((200, 200), PIL.Image.ANTIALIAS)
                    timage.save(
                        os.path.join(
                            '%s-zoomed' % os.path.basename(self.series[series_index].file_path),
                            'slide_thumbnail.png'
                        )
                    )
                    # Create large thumbnail
                    timage = image.copy()
                    timage.thumbnail((1000, 1000), PIL.Image.ANTIALIAS)
                    timage.save(
                        os.path.join(
                            '%s-zoomed' % os.path.basename(self.series[series_index].file_path),
                            'slide_large_thumbnail.png'
                        )
                    )
        return failed_message

    def load_series(self, i_xyzct=(-1,-1,-1,-1,-1), n_xyzct=(-1,-1,-1,-1,-1)):
        """
        Load series of images from study directory. Once it determines necessary information such as version information
        and a list of all TIFF files in the study directory, the processor then loads each series image,
        sets image dimension properties and generates exploded view assets. No exceptions are explicity raised but
        catches multiple IOError and StopException exceptions and continues to process rest of images.

        Adapted from https://gist.github.com/jzaremba/f7adde4dced95cbdbf9b

        TODO: Implement for multi-plane (z stack) studies. Currently, only supports z=1

        Args:

        Returns:

        Raises:

        """
        tif_file_regex = re.compile('\S_Cycle.*Ch[12]_0+\d+.*tif')
        failed_messages = []

        for image_directory, folders, files in os.walk(self.unzipped_file_dir):
            # Assemble a list of all .tif files
            tif_files = [f for f in files if f.endswith('.ome.tif')]

            # Ensure list of tif files is legitimate
            if len(tif_files) < 3:
                err_msg = 'Less than 3 files in {}'.format(image_directory)
                self.LOGGER.error(err_msg)
                failed_messages.append(err_msg)
                continue
            self.LOGGER.info('Discovered directory of tiff files in %s' % image_directory)

            # Determine filename of config XML file
            xml_files = [f for f in files if f.endswith('.xml')]
            if len(xml_files) > 1:
                err_msg = 'Too many XML files in directory {}'.format(image_directory)
                self.LOGGER.error(err_msg)
                failed_messages.append(err_msg)
                continue

            try:
                match = re.search(tif_file_regex, tif_files[0])
                basename = tif_files[0][:match.start() + 1]
                xml_file_name = basename + '.xml'
            except Exception:
                xml_file_name = xml_files[0]
                err_msg = 'No XML file in {} matches standard tif ' \
                          'file naming convention. Using {} instead.'.format(image_directory, xml_file_name)
                self.LOGGER.warn(err_msg)
                failed_messages.append(err_msg)

            # Lock the xml file, if it can't be locked just continue
            # This will hold the lock if it exits early before conversion,
            # but it will be released again when the process ends
            lockfile = open(os.path.join(image_directory, xml_file_name))
            try:
                fcntl.flock(lockfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except IOError:
                err_msg = 'XML header file {} in directory {} ' \
                          'cannot be locked for reading and processing.'.format(xml_file_name, image_directory)
                self.LOGGER.error(err_msg)
                failed_messages.append(err_msg)
                continue

            # Parse the xml name, Prairie version, and cfg filename
            try:
                prairie_version, cfg_filename, protocol_elements, cycles = self.parse_prairieview_xml(
                    os.path.join(image_directory, xml_file_name)
                )
            except Exception:
                err_msg = 'Exception when parsing xml {} ' \
                          'to get Prairie version and config ' \
                          'cfg file with exception info: {}'.format(xml_file_name, sys.exc_info())
                self.LOGGER.error(err_msg)
                failed_messages.append(err_msg)
                continue

            # Create a generator of all the Sequences in the XML
            sequences = (
                elem for _, elem in
                ElementTree.iterparse(os.path.join(image_directory, xml_file_name))
                if elem.tag == 'Sequence' and elem.get('type') != 'TSeries Voltage Output Experiment'
            )

            # Iterate over protocols and cycles
            # There should be one OME-TIFF file per cycle-protocol.
            tiffs_to_save, failed_message = self.collect_tiffs_to_save(
                cycles,
                protocol_elements,
                sequences,
                prairie_version,
                os.path.join(image_directory, xml_file_name)
            )

            # Check tiffs collected to decode and log errors
            try:
                self.check_tiffs_to_save(tiffs_to_save, sequences, failed_message)
            except Exception:
                err_msg = 'Exception when checking tiffs in directory {}: {}'.format(image_directory, sys.exc_info())
                self.LOGGER.error(err_msg)
                failed_messages.append(err_msg)
                continue

            failed_message = self.save_outputs(
                tiffs_to_save,
                os.path.join(image_directory, xml_file_name),
                prairie_version,
                i_xyzct=i_xyzct,
                n_xyzct=n_xyzct
            )

            # Log any failures encountered during processor workflow
            if failed_message:
                self.LOGGER.error(traceback.print_tb(failed_message[2]))
                failed_messages.append(failed_message)

            # Release the file lock
            lockfile.close()
        return failed_messages

    def task(self):
        """ Run processor task. Generates output images (OME-TIFF and exploded 'view' assets) and publishes outputs. """

        # self._load_image()
        self.LOGGER.info('Got inputs {}'.format(self.inputs))

        # Get sub_region index
        try:
            sub_region = self.inputs['sub_region_file']
            sub_region_regex = r'sub_' \
                               r'x_([0-9]+)_([0-9]+)_' \
                               r'y_([0-9]+)_([0-9]+)_' \
                               r'z_([0-9]+)_([0-9]+)_' \
                               r'c_([0-9]+)_([0-9]+)_' \
                               r't_([0-9]+)_([0-9]+).txt'
            i_x = int(re.match(re.compile(sub_region_regex), sub_region).groups()[0])
            n_x = int(re.match(re.compile(sub_region_regex), sub_region).groups()[1])
            i_y = int(re.match(re.compile(sub_region_regex), sub_region).groups()[2])
            n_y = int(re.match(re.compile(sub_region_regex), sub_region).groups()[3])
            i_z = int(re.match(re.compile(sub_region_regex), sub_region).groups()[4])
            n_z = int(re.match(re.compile(sub_region_regex), sub_region).groups()[5])
            i_c = int(re.match(re.compile(sub_region_regex), sub_region).groups()[6])
            n_c = int(re.match(re.compile(sub_region_regex), sub_region).groups()[7])
            i_t = int(re.match(re.compile(sub_region_regex), sub_region).groups()[8])
            n_t = int(re.match(re.compile(sub_region_regex), sub_region).groups()[9])
        except (KeyError, IndexError, AttributeError):
            i_x = -1
            n_x = -1
            i_y = -1
            n_y = -1
            i_z = -1
            n_z = -1
            i_c = -1
            n_c = -1
            i_t = -1
            n_t = -1

        # Unzip file asset
        os.system('tar -zxvf %s' % self.file)
        assert os.path.isdir('%s' % os.path.basename(self.file).replace('.brukertiff.gz', ''))
        # Clean zip asset
        os.system('rm %s' % self.file)
        os.system('rm -r */*MIP*')
        os.system('rm */RoiSet.zip')
        os.system('rm */*.xlsx')

        self.unzipped_file_dir = os.path.basename(self.file).replace('.brukertiff.gz', '')

        # Load and save view images
        failed_messages = self.load_series(i_xyzct=(i_x, i_y, i_z, i_c, i_t), n_xyzct=(n_x, n_y, n_z, n_c, n_t))
        if failed_messages:
            for failed_message in failed_messages:
                self.LOGGER.error('The following error was raised during execution of '
                                  'brukertiff-processor: {}'.format(failed_message))

        for series_index in range(len(self.series)):
            # Get output OMETIFFFile
            output_file = self.series[series_index]
            self.LOGGER.info('Current list of files include: {}'.format(os.listdir('.')))

            # Save dimensions object as JSON in view/ directory (for now)
            with open(os.path.join('%s-zoomed' % os.path.basename(output_file.file_path), 'dimensions.json'),
                      'w') as fp:
                json.dump(output_file.img_dimensions, fp)

            # Create create-asset JSON object file called view_asset_info.json
            self.upload_key = os.path.join(
                self.settings.storage_directory,
                os.path.basename(output_file.file_path) + '-zoomed'
            )
            with open('view_asset_info.json', 'w') as fp:
                json.dump(output_file.get_view_asset_dict(
                    self.settings.storage_bucket,
                    self.upload_key
                ), fp)

            # Generate properties metadata.json metadata
            metadata = []
            img_dimensions = output_file.img_dimensions
            for dim in range(output_file.num_dimensions):
                for property_key_suffix in ["assignment", "length", "resolution", "units"]:
                    # Initialize property
                    property = {}

                    # Set property key and value
                    property_key = 'dimensions.%i.%s' % (dim, property_key_suffix)
                    property_value = str(img_dimensions['dimensions'][dim][property_key_suffix])

                    # Create property instance
                    property["key"] = property_key
                    property["value"] = property_value
                    property["dataType"] = "String"
                    property["category"] = "Blackfynn"
                    property["fixed"] = False
                    property["hidden"] = False
                    metadata.append(property)
            with open('metadata.json', 'w') as fp:
                json.dump(metadata, fp)
