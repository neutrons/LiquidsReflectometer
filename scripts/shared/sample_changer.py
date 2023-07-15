import sys
import time
from epics import PV

MOVE_TIMEOUT = 600
ALIGN_TIMEOUT = 1200

# Set to True only for debugging 
RETURN_ON_FAIL = True

# Data types: 0=data, 1=direct beam (full), 2=direct beam (0-att), 3=ignore
DATATYPE = PV('BL4B:CS:ExpPl:DataType')
AUTOALIGN = PV('BL4B:CS:AutoAlign:Run')
AUTOALIGN_STATUS = PV('BL4B:CS:AutoAlign:Stat')

# 1 is running, 0 is idle
RUNNING = PV('BL4B:CS:Running:Stat')
BL4B_MOT_PREFIX = 'BL4B:Mot:'

POSITIONS = [dict(ysc=50,  hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=108, hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=165, hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=222, hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=279, hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=336, hs=9, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=50,  hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=108, hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=165, hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=222, hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=279, hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             dict(ysc=336, hs=66, ths=-3.7, chis=0, zs=199.2, xs=0, ys=6.002),
             ]

SCAN = [[6000, {'BL4B:Chop:Gbl:WavelengthReq': 15,     's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0.6, 'tthd': 1.2, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [6000, {'BL4B:Chop:Gbl:WavelengthReq': 12.386, 's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0.6, 'tthd': 1.2, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [6000, {'BL4B:Chop:Gbl:WavelengthReq': 9.74,   's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0.6, 'tthd': 1.2, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [6000, {'BL4B:Chop:Gbl:WavelengthReq': 7.043,  's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0.6, 'tthd': 1.2, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [18000, {'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 0.39, 'si:Y:Gap': 0.25, 's3:Y:Gap': 10, 'ths': 0.6, 'tthd': 1.2, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [18000, {'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 0.769, 'si:Y:Gap': 0.493, 's3:Y:Gap': 10, 'ths': 1.183, 'tthd': 2.366, 's1:X:Gap': 20, 'si:X:Gap': 20}],
        [72000, {'BL4B:Chop:Gbl:WavelengthReq': 4.25,  's1:Y:Gap': 1.523, 'si:Y:Gap': 0.976, 's3:Y:Gap': 20, 'ths': 2.343, 'tthd': 4.686, 's1:X:Gap': 20, 'si:X:Gap': 20}],
       ]


def move_motors(positions):
    check_list = []
    print("Moving:")
    for motor, position in enumerate(positions):
        print("  %s -> %s" % (motor, position))
        if motor.startswith("BL4B"):
            _motor = motor
        else:
            _motor = BL4B_MOT_PREFIX + motor
        _pv = PV(_motor).put(position)
        check_list.append(PV(_motor + '.Status'))

    ready = False
    t0 = time.time()
    while not ready:
        time.sleep(2)
        print('  ..checking')
        for _pv in check_list:
            # Check motor status
            ready = ready and _pv.get() == 0
        if time.time() - t0 > MOVE_TIMEOUT:
            print("Timed out....")
            return RETURN_ON_FAIL

    print('Ready')
    return True


def align_sample():
    print("Starting automated alignment")
    t0 = time.time()
    AUTOALIGN.put(1)
    while AUTOALIGN.get() is not 1:
        if time.time() - t0 > ALIGN_TIMEOUT:
            print("Timed out...")
            return RETURN_ON_FAIL
        print("  ..checking")
        time.sleep(10)
    print("Sample is align")
    return True


def run_scan(name, scan):
    group_id = PV('BL4B:CS:RunControl:LastRunNumber').get() + 1

    PV("BL4B:CS:Autoreduce:Sequence:Total").put(len(scan))
    PV("BL4B:CS:Autoreduce:Sequence:Id").put(group_id)
    sequence_num = PV("BL4B:CS:Autoreduce:Sequence:Num")
    title = PV('BL4B:CS:Autoreduce:BaseTitle')

    for i in range(1, len(scan)+1):
        sequence_num.put(i)
        title.put("%s-%s-%s" % (name, group_id, i))

        if move_motors(scan[i][1]):
            # Acquire neutrons!
            PV('BL4B:CS:RunControl:Start').put(1)
            neutrons=PV('BL4B:Det:Neutrons')
            while neutrons.get() < scan[i][0]:
                time.sleep(5)
            PV('BL4B:CS:RunControl:Stop').put(1)

    return True


def process(positions):
    DATATYPE.put(0)

    for sample in positions:
        # Move the sample changer
        assert(sample['pos']>0 and sample['pos']<13)
        if move_motors(sample['pos']+1):
            # Align the sample
            if align_sample():
                # Run the CSV
                run_scan(sample['name'], SCAN)


if __name__ == "__main__":
    positions = [dict(name='quartz', pos=1),
                 dict(name='SiO2', pos=2)
                 ]
    process(positions)