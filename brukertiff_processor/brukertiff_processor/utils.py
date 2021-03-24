import fileinput
import sys
from datetime import datetime as dt
from distutils.version import LooseVersion
from xml.etree import ElementTree


def timestamp():
    return dt.strftime(dt.now(), '[%Y-%m-%d %H:%M:%S]')


def get_element_size_um(xml_file_path, prairie_version):
    """
    Get iamge resolution and units for any plane in Bruker PrairieView image set.
    TODO: Implement setting multi-plane resolution (in z dimension)

    Args:
        xml_file_path: Absolute path to XML-based config file for image set.
        prairie_version: str of version of Prairie config file.

    Returns:
        Tuple containing resolution for (z, x, y)

    Raises:
        None
    """

    # Initialize dimension resolution
    x = 1  # um
    y = 1  # um

    if prairie_version >= LooseVersion('5.2'):
        for _, elem in ElementTree.iterparse(xml_file_path):
            if elem.get('key') == 'micronsPerPixel':
                for value in elem.findall('IndexedValue'):
                    if value.get('index') == "XAxis":
                        x = float(value.get('value'))
                    elif value.get('index') == "YAxis":
                        y = float(value.get('value'))
                return 1, y, x
    else:
        for _, elem in ElementTree.iterparse(xml_file_path):
            if elem.tag == 'PVStateShard':
                for key in elem.findall('Key'):
                    if key.get('key') == 'micronsPerPixel_XAxis':
                        x = float(key.get('value'))
                    elif key.get('key') == 'micronsPerPixel_YAxis':
                        y = float(key.get('value'))
                return 1, y, x

    # If not found, assume 1 um x 1 um x 1 um
    return 1, y, x


def get_prairieview_version(xml_file_path):
    """
    Returns the version number of the PrairieView Bruker imaging dataset by reading the XML config file.
    Adapted from https://gist.github.com/jzaremba/f7adde4dced95cbdbf9b.

    Args:
        xml_file_path: XML configuration file.

    Returns:
        Returns a str containing the version number.

    Raises:
        None
    """

    for _, elem in ElementTree.iterparse(xml_file_path, events=("start",)):
        if elem.tag == 'PVScan':
            return LooseVersion(elem.get('version'))


def parse_cfg_file(xml_file_path):
    """
    Parses Bruker ParirieView config file and determines the number of experiments and iterations.
    Adapted from https://gist.github.com/jzaremba/f7adde4dced95cbdbf9b

    Args:
        xml_file_path: Absolute file path to config file for PrairieView file (should be XML).

    Returns:
        Returns a tuple containing all time or plane series along with number of iterations.

    Raises:
        None
    """

    try:
        cfg_tree = ElementTree.parse(xml_file_path)
        # Shouldn't need to except the TypeError, but something is failing
        # within the ElementTree parser
    except (ElementTree.ParseError, TypeError):
        reformat_prairie_cfg(xml_file_path)
        cfg_tree = ElementTree.parse(xml_file_path)

    tSeries_element = cfg_tree.find('TSeries')
    try:
        n_iterations = int(tSeries_element.get('iterations'))
    except TypeError:
        n_iterations = None

    elements = []
    for elem in tSeries_element:
        if elem.tag in ['PVTSeriesElementSequence', 'PVTSeriesElementZSeries']:
            elements.append((elem.tag, int(elem.get('repetitions'))))

    return elements, n_iterations


def reformat_prairie_cfg(xml_file_path):
    """
    Reformat configuration file from older PrairieView files (*.cfg) by replacing all newlines with spaces.
    Adapted from https://gist.github.com/jzaremba/f7adde4dced95cbdbf9b.

    Args:
        xml_file_path: Absolute file path to config file for PrairieView file.

    Returns:
        Returns the same config file with replaced characters

    Raises:
        None
    """

    for line in fileinput.input(xml_file_path, inplace=1):
        if '&#x1;' in line:
            line = line.replace('&#x1;', ' ')
        sys.stdout.write(line)
