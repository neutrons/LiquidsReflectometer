"""
    Runs diagnostics to plan a set of composite direct beams with Cd.
"""
import argparse
import csv
import json

import instrument

#SCAN = [
#    [100, {'BL4B:Chop:Gbl:SpeedReq': 30, 'BL4B:Chop:Gbl:WavelengthReq': 6, 's1:Y:Gap': 10.528, 'si:Y:Gap': 1.128, 's3:Y:Gap': 30, 'ths': 0., 'tthd': 0, 's1:X:Gap': 20, 'si:X:Gap': 20, 
#           'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':0, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}, (4,4), 'C-DB 30Hz angle 3'],
#        ]

ACTUATORS = {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1}

# Order to put the actuators back in
TRYOUT_ORDER = [
    {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':0, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':0, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':0, 'BL4B:Actuator:100M':0, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':0, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':0, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':0, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':0, 'BL4B:Actuator:200M':0, 'BL4B:Actuator:400M':1},
    {'BL4B:Actuator:50M':1, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':0},
    {'BL4B:Actuator:50M':0, 'BL4B:Actuator:100M':1, 'BL4B:Actuator:200M':1, 'BL4B:Actuator:400M':0},
]

def read_scan_configuration(scan_file: str):
    """
    Reads the scan configuration from a file.
    
    :param scan_file: The scan configuration file.
    :return: The scan configuration.
    """
    data = []
    with open(scan_file, mode='r', newline='') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.append(row)
    return data


def remove_chars_before_ref(s: str, ref_char: str) -> str:
    """
    Removes characters in a string in front of a reference character.
    
    :param s: The input string.
    :param ref_char: The reference character.
    :return: The substring starting from the reference character.
    """
    ref_index = s.find(ref_char)
    if (ref_index != -1):
        return s[ref_index:]
    return s


def process_configuration_list(config: list):
    """
    Processes a configuration list and return a dictionary of PV values
    that can be used to plan the composite direct beams.
    
    :param config: The list of configurations.
    :return: The list of PV values.
    """
    pv_list = config[0]
    scan_list = []

    for i, item in enumerate(config[1:]):
        pv_dict = {}
        for j, pv in enumerate(pv_list):
            if 'BL4B' not in pv or 'Autoreduce' in pv:
                continue
            # Remove annotations in front of BL4B in each PV name
            _pv_name = remove_chars_before_ref(pv, 'BL4B')
            _pv_name = _pv_name.replace(instrument.BL4B_MOT_PREFIX, '')
            pv_dict[_pv_name] = float(item[j])
            # Set the angles to zero
            if 'ths' in _pv_name or 'tthd' in _pv_name:
                pv_dict[_pv_name] = 0

        scan_list.append(pv_dict)
    return scan_list


def find_best_n_for_slicing(scan: dict, rate_cutoff: int, n_cutoff: int):
    """
    Finds the best N for slicing the slits.
    
    :param scan: The scan configuration.
    :param rate_cutoff: The rate cutoff for the scan in neutrons/second.
    :param n_cutoff: The cutoff for the number of slit slices.
    """
    lr = instrument.LiquidsReflectometer()
    si_x_gap = float(scan['si:X:Gap'])
    s1_x_gap = float(scan['s1:X:Gap'])
    try_si_x_gap = si_x_gap
    try_s1_x_gap = s1_x_gap

    best_n = -1
    for n in range(1, n_cutoff):
        try_si_x_gap = si_x_gap / n
        try_s1_x_gap = s1_x_gap / n
        lr.move({'si:X:Gap': try_si_x_gap, 's1:X:Gap': try_s1_x_gap})

        # Measure neutron rate
        rate = lr.measure_rate(5)
        print("  Measured: %s; rate=%g" % (n, rate))

        # Adjust N if rate is above cutoff

        # Break if rate is below cutoff
        if rate < rate_cutoff:
            best_n = n
            break
    return best_n


def plan_composite_direct_beams(scan_file: str, rate_cutoff=10000, n_cutoff=5,
                                charge:float = 100, title: str = ''):
    """
    Plans a set of composite direct beams with Cd.
    
    :param scan_file: The scan configuration file.
    :param rate_cutoff: The rate cutoff for the scan in neutrons/second.
    :param n_cutoff: The cutoff for the number of slit slices.
    :param title: The title of the scan.
    """
    data = read_scan_configuration(scan_file)

    scan_list = process_configuration_list(data)
    final_list = []

    lr = instrument.LiquidsReflectometer()
 
    for i, scan in enumerate(scan_list):
        print("Scan configuration: ", i)
        # Move to position
        lr.move(scan)
        # Move out actuators
        lr.move(ACTUATORS)

        # Measure neutron rate and adjust N
        best_n = find_best_n_for_slicing(scan, rate_cutoff, n_cutoff)

        # If we didn't find a good N, add Cd
        tryout = ACTUATORS
        actuator_pos = -1
        if best_n < 0:
            move_actuators = False
            si_x_gap = float(scan['si:X:Gap']) / n_cutoff
            s1_x_gap = float(scan['s1:X:Gap']) / n_cutoff
            lr.move({'si:X:Gap': si_x_gap, 's1:X:Gap': s1_x_gap})

            for j, tryout in enumerate(TRYOUT_ORDER):
                lr.move(tryout)
                rate = lr.measure_rate(5)
                print("  Measured: %s; %s; rate=%g" % (n_cutoff, j, rate))
                if rate < rate_cutoff:
                    move_actuators = True
                    actuator_pos = j
                    break

            if actuator_pos < 0:
                print("Config {i}: COULD NOT FIND LOW-RATE CONFIGURATION".format(i=i))
                print(r'Config {i}: N={best_n} Cd={tryout}'.format(i=i, best_n=best_n, tryout=tryout))
            else:
                scan.update(tryout)
                final_list.append([charge, scan, (n_cutoff, n_cutoff), "%s N=%s Cd=%s" % (title, n_cutoff, actuator_pos)])
        else:
            final_list.append([charge, scan, (best_n, best_n), "%s N=%s Cd=None" % (title, best_n)])

    # Write list as json file
    with open('composite_direct_beams_%s.json' % title, 'w') as file:
        json.dump(final_list, file)

    print(final_list)
    return final_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plan composite direct beams with Cd.")
    parser.add_argument("--scan", help="The scan configuration file to use.")
    parser.add_argument("--title", help="The title of the scan [like the medium].", default='')
    parser.add_argument("--rate_max", help="Max neutrons/sec", default=10000, type=float)
    parser.add_argument("--n_max", help="Max N", default=5, type=int)
    parser.add_argument("--charge", help="The charge to use for the scan.",
                        default=100, type=float)
    args = parser.parse_args()

    print("Planning composite direct beams with Cd for the reflectometer.")
    print("Scan configuration file: ", args.scan)
    print("Title: ", args.title)
    plan_composite_direct_beams(
        args.scan, charge=args.charge, rate_cutoff=args.rate_max,
        n_cutoff=args.n_max, title=args.title
    )
    print("Planned composite direct beams with Cd.")
    print("Done.")

