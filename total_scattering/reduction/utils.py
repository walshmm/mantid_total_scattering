import itertools
import numpy as np
from mantid.simpleapi import CreateEmptyTableWorkspace


def expand_ints(s):
    ''' Function to expand string of ints with dashes
        Ex. "1-3, 8-9, 12" -> [1,2,3,8,9,12]

    :param s: string with the ints we want to expand
    :type s: str

    :return: List of integers from expansion
    :rtype: list[int]
    '''
    spans = (el.partition('-')[::2] for el in s.split(','))
    ranges = (range(int(s), int(e) + 1 if e else int(s) + 1)
              for s, e in spans)
    all_nums = itertools.chain.from_iterable(ranges)
    return list(all_nums)


def compress_ints(int_list):
    ''' Function to compress list of ints with dashes
    Ex. [1,2,3,8,9,12] -> 1-3, 8-9, 12

    :param int_list: List of integers to compress
    :type int_list: list[int]

    :return: String for the "compressed" integers
    :rtype: str
    '''
    seq = []
    final = []
    last = 0

    for index, val in enumerate(int_list):

        if last + 1 == val or index == 0:
            seq.append(val)
            last = val
        else:
            if len(seq) > 1:
                final.append(str(seq[0]) + '-' + str(seq[len(seq) - 1]))
            else:
                final.append(str(seq[0]))
            seq = []
            seq.append(val)
            last = val

        if index == len(int_list) - 1:
            if len(seq) > 1:
                final.append(str(seq[0]) + '-' + str(seq[len(seq) - 1]))
            else:
                final.append(str(seq[0]))

    final_str = ', '.join(map(str, final))
    return final_str


def extract_key_match_from_dict(keys, dictionary):
    """ Convienence function for extraction of one key from dictionary

    :param keys: Keys to check against dictionary
    :type keys: list
    :param dictionary: Dictionary to check
    "type dictionary: dict

    :return: The exctracted value
    :rtype: any
    """
    out = find_key_match_in_dict(keys, dictionary)
    if out:
        return out
    else:
        e = "No matching key found. Valid keys are {}".format(keys)
        raise Exception(e)


def find_key_match_in_dict(keys, dictionary):
    """ Check if one and only one of the keys is in the dictionary
    and return its value

    :param key: Keys we will check for in dictionary
    :type key: str
    :param dictionary: Dictionary to check
    :type dictionary: dict

    :return: Either the value in dictionary for the key or None if not found
    :rtype: value in dict or None
    """
    # Get the boolean for each key if it exists in the dictionary
    keys_exist_in_dict = map(lambda key: key in dictionary, keys)

    # If only one exists, return the match, else raise exception
    if one_and_only_one(keys_exist_in_dict):
        for key in keys:
            if key in dictionary:
                return dictionary[key]

    # None of the keys in the dictionary, return None
    return None


def generate_cropping_table(qmin, qmax):
    ''' Generate a Table workspace that can be used to crop
    another workspace in reciprocal space (ie MomentumTransfer)

    :param qmin: list of Qmin values to crop each spectra by
    :type qmin: str (example: '0.2,0.4,1.0')
    :param qmax: list of Qmax values to crop each spectra by
    :type qmax: str (example: '10.5,12.0,40.0')

    :return: Cropping table with columns of ("SpectraList","Xmin","Xmax")
    :rtype: TableWorkspace
    '''
    mask_info = CreateEmptyTableWorkspace()
    mask_info.addColumn("str", "SpectraList")
    mask_info.addColumn("double", "XMin")
    mask_info.addColumn("double", "XMax")
    for (i, value) in enumerate(qmin):
        mask_info.addRow([str(i), 0.0, value])
    for (i, value) in enumerate(qmax):
        mask_info.addRow([str(i), value, 100.0])

    return mask_info


def get_each_spectra_xmin_xmax(wksp):
    ''' Get Xmin and Xmax lists for Workspace, excluding
    values of inf and NaN

    :param wksp: Workspace to extract the xmin and xmax values from
    :type qmin: Mantid Workspace

    :return: Lists for XMin and XMax values:
    :rtype: (list, list) == (xmin, xmax)
    '''
    xmin = list()
    xmax = list()
    numSpectra = wksp.getNumberHistograms()
    for i in range(numSpectra):
        x = wksp.readX(i)
        xmin.append(np.nanmin(x[x != -np.inf]))
        xmax.append(np.nanmax(x[x != np.inf]))
    return xmin, xmax



def get_normalization(config):
    """ Extract the normalization section from JSON input

    :param config: JSON input for reduction
    :type config: dict

    :return: The exctracted value for normalization in the input
    :rtype: any
    """
    keys = ["Normalization", "Normalisation", "Vanadium"]
    out = extract_key_match_from_dict(keys, config)
    return out


def get_sample(config):
    """ Extract the sample section from JSON input

    :param config: JSON input for reduction
    :type config: dict

    :return: The exctracted value for sample in the input
    :rtype: any
    """
    keys = ["Sample"]
    out = extract_key_match_from_dict(keys, config)
    return out


def one_and_only_one(iterable):
    """Determine if iterable (ie list) has one and only one `True` value

    :param iterable: The iterable to check
    :type iterable: list

    :return: If there is one and only one True
    :rtype: bool
    """
    try:
        iterator = iter(iterable)
        has_true = any(iterator)
        has_another_true = any(iterator)
        return has_true and not has_another_true
    except Exception as e:
        print(e)
        raise
