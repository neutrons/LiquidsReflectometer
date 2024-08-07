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
StartDiag = PV("BL4B:Det:N1:Start")
StopDiag = PV("BL4B:Det:N1:Stop")
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
RETURN_ON_FAIL = False


class LiquidsReflectometer:

    def __init__(self, is_virtual: bool = False):
        self.is_virtual = is_virtual
        if not is_virtual and IS_VIRTUAL:
            print("Warning: EPICS module not found. Running in virtual mode.")
            self.is_virtual = True

        # Stop any existing runs
        StopRun.put(1)
        StopDiag.put(1)
        self.acquiring = False

        # Virtual values
        self.virtual_counts = 0
        self.virtual_timer = 0
    
    def initialize_series(self, length: int = 1, title='Composite DB'):
        group_id = PV('BL4B:CS:RunControl:LastRunNumber').get() + 1
        PV("BL4B:CS:Autoreduce:Sequence:Num").put(1)
        PV("BL4B:CS:Autoreduce:Sequence:Id").put(group_id)
        PV('BL4B:CS:Autoreduce:BaseTitle').put(title)
        PV("BL4B:CS:Autoreduce:Sequence:Total").put(length)
        PV("BL4B:CS:Autoreduce:DataType").put(3)
        #PV('BL4B:CS:ExpPl:DataType').put(5)

    def increment_sequence(self, title='C-DB'):
        PV('BL4B:CS:Autoreduce:BaseTitle').put(title)
        sequence_num = PV("BL4B:CS:Autoreduce:Sequence:Num")
        sequence_num.put(sequence_num.get() + 1)

    def move(self, positions):
        check_list = []
        sleep_for_actuators = False
        print("Moving:")
        for i, (motor, position) in enumerate(positions.items()):
            if 'ths' in motor:
                continue
            print("  %s -> %s" % (motor, position))
            if motor.startswith("BL4B"):
                _motor = motor
            else:
                _motor = BL4B_MOT_PREFIX + motor
            _pv = PV(_motor).put(position)
            if 'Wavelength' in _motor:
                check_list.append(PV('BL4B:Chop:Gbl:Busy:Stat'))
            elif 'SpeedReq' in _motor:
                pass
            elif 'Actuator' in _motor:
                sleep_for_actuators = True
            else:
                check_list.append(PV(_motor + ':Status'))

        # Sleep for actuators
        if sleep_for_actuators:
            time.sleep(5)

        ready = IS_VIRTUAL
        t0 = time.time()
        while not ready:
            time.sleep(0.5)
            #print('  ... checking')
            ready = True
            for _pv in check_list:
                # Check motor status
                status = _pv.get()
                #print(type(status), status)
                ready = ready and _pv.get() == 0
            if time.time() - t0 > MOVE_TIMEOUT:
                print("Timed out ...")
                return RETURN_ON_FAIL

        print('Ready')
        return True
    
    def start_or_resume(self, counts: int = 0, seconds: int = 0, charge: float = 200):
        """
            Start or resume a run depending on the current state.
        """
        print("Acquire [current state: %s] %g %g %g" % (self.acquiring, counts, seconds, charge))
        if self.acquiring:
            PauseRun.put(0)
            # Get the current charge
            _c = C.get()
        else:
            StartRun.put(1)
            self.acquiring = True
            # Make sure we start at zero charge. Doing this now
            # prevents a race condition.
            _c = 0

        # If we are virtual, update the virtual counts and timer
        if self.is_virtual:
            self.virtual_counts = counts if counts > 0 else 1000
            self.virtual_timer = seconds if seconds > 0 else 1000
            return



        time.sleep(1)
        # Wait for the neutron count to reach the desired value
        if counts > 0:
            while neutrons.get() < counts:
                time.sleep(0.1)
        # Wait for the desired number of seconds
        elif seconds > 0:
            time.sleep(seconds)
        elif charge > 0:
            charge_to_reach = charge + _c
            print("Starting charge: %s -> %s" % (_c, charge_to_reach))
            while C.get() < charge_to_reach:
                #TODO: may want to consider the following to ensure that the DAS is still
                # running. If it's not we should just exit. It may mean someone clicked StopAll
                #is_running = PV("BL4B:CS:RunControl:Running").get()
                #if is_running == 0:
                #    sys.exit(0)
                _c_check = C.get()
                #print("    q=%s" % _c_check)
                time.sleep(0.1)
    
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
    
    def get_deadtime(self):
        print("Not yet implemented")

    
    def get_rate(self):
        """
            Get the current count rate.
        """
        if self.is_virtual:
            return self.virtual_counts / self.virtual_timer

        total_neutrons = neutrons.get()
        total_time = timer.get()
        return total_neutrons / total_time

    def measure_rate(self, time_interval: int = 10):
        """
        Measure the rate for a given time.
        """
        StartDiag.put(1)
        time.sleep(time_interval)
        StopDiag.put(1)

        total_neutrons = neutrons.get()
        total_time = timer.get()
        return total_neutrons / total_time
