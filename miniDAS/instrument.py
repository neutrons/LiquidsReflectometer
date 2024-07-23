import time

# Import EPICS module but allow for a virtual environment
IS_VIRTUAL = False
try:
    from epics import PV
except ImportError:
    IS_VIRTUAL = True
    class PV:
        def __init__(self, pvname):
            self.name = pvname
            self.value = -1
        def put(self, value):
            self.value = value
        def get(self):
            return self.value

# Define channel access variables
S1X = PV('BL4B:Mot:s1:X:Gap')
S1Y = PV('BL4B:Mot:s1:Y:Gap')
S1Xrbv = PV('BL4B:Mot:s1:X:Gap:Readback')
S1Yrbv = PV('BL4B:Mot:s1:Y:Gap:Readback')
SiX = PV('BL4B:Mot:si:X:Gap')
SiY = PV('BL4B:Mot:si:Y:Gap')
SiXrbv = PV('BL4B:Mot:si:X:Gap:Readback')
SiYrbv = PV('BL4B:Mot:si:Y:Gap:Readback')

S1XCEN = PV('BL4B:Mot:s1:X:Center')
S1YCEN = PV('BL4B:Mot:s1:Y:Center')
S1XCENrbv = PV('BL4B:Mot:s1:X:Center:Readback')
S1YCENrbv = PV('BL4B:Mot:s1:Y:Center:Readback')

SiXCEN = PV('BL4B:Mot:si:X:Center')
SiYCEN = PV('BL4B:Mot:si:Y:Center')
SiXCENrbv = PV('BL4B:Mot:si:X:Center:Readback')
SiYCENrbv = PV('BL4B:Mot:si:Y:Center:Readback')

C = PV('BL4B:Det:PCharge')
neutrons = PV('BL4B:Det:Neutrons')
timer = PV('BL4B:CS:RunControl:RunTimer')
StartRun = PV('BL4B:CS:RunControl:Start')
StopRun = PV('BL4B:CS:RunControl:Stop')
PauseRun = PV('BL4B:CS:RunControl:Pause')
Rate = PV('BL4B:Det:N1:Det1:EventRate_RBV')
Name = PV('BL4B:CS:Autoreduce:BaseTitle')

Lcen = PV('BL4B:Det:TH:BL:Lambda')
Lset = PV('BL4B:Chop:Gbl:WavelengthReq')
ChopStat = PV('BL4B:Chop:Gbl:Busy:Stat')

BL4B_MOT_PREFIX = 'BL4B:Mot:'

MOVE_TIMEOUT = 600

# Set to True only for debugging 
RETURN_ON_FAIL = True


class LiquidsReflectometer:

    def __init__(self, is_virtual: bool = False):
        self.is_virtual = is_virtual
        if not is_virtual and IS_VIRTUAL:
            print("Warning: EPICS module not found. Running in virtual mode.")
            self.is_virtual = True

        # Stop any existing runs
        StopRun.put(1)
        self.acquiring = False

        # Virtual values
        self.virtual_counts = 0
        self.virtual_timer = 0
    
    def initialize_series(self, seq: int = 1, length: int = 1, title='Composite DB'):
        group_id = PV('BL4B:CS:RunControl:LastRunNumber').get() + 1
        PV("BL4B:CS:Autoreduce:Sequence:Total").put(length)
        PV("BL4B:CS:Autoreduce:Sequence:Id").put(group_id)
        PV('BL4B:CS:Autoreduce:BaseTitle').put(title)
        PV('BL4B:CS:ExpPl:DataType').put(5)

    def increment_sequence(self):
        sequence_num = PV("BL4B:CS:Autoreduce:Sequence:Num")
        sequence_num.put(sequence_num.get() + 1)

    def move(self, positions):
        check_list = []
        print("Moving:")
        for i, (motor, position) in enumerate(positions.items()):
            print("  %s -> %s" % (motor, position))
            if motor.startswith("BL4B"):
                _motor = motor
            else:
                _motor = BL4B_MOT_PREFIX + motor
            _pv = PV(_motor).put(position)
            check_list.append(PV(_motor + '.Status'))

        ready = IS_VIRTUAL
        t0 = time.time()
        while not ready:
            time.sleep(2)
            print('  ... checking')
            for _pv in check_list:
                # Check motor status
                ready = ready and _pv.get() == 0
            if time.time() - t0 > MOVE_TIMEOUT:
                print("Timed out ...")
                return RETURN_ON_FAIL

        print('Ready')
        return True
    
    def start_or_resume(self, counts: int = 0, seconds: int = 0):
        """
            Start or resume a run depending on the current state.
        """
        if self.acquiring:
            PauseRun.put(0)
        else:
            StartRun.put(1)
            self.acquiring = True

        # If we are virtual, update the virtual counts and timer
        if self.is_virtual:
            self.virtual_counts = counts if counts > 0 else 1000
            self.virtual_timer = seconds if seconds > 0 else 1000
            return

        # Wait for the neutron count to reach the desired value
        if counts > 0:
            while neutrons.get() < counts:
                time.sleep(2)
        
        # Wait for the desired number of seconds
        elif seconds > 0:
            time.sleep(seconds)
    
    def pause(self):
        """
            Pause the current run.
        """
        PauseRun.put(1)

    def stop(self):
        """
            Stop the current run.
        """
        StopRun.put(1)
        self.acquiring = False
    
    def get_rate(self):
        """
            Get the current count rate.
        """
        if self.is_virtual:
            return self.virtual_counts / self.virtual_timer

        total_neutrons = neutrons.get()
        total_time = timer.get()
        return total_neutrons / total_time
