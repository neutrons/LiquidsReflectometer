import argparse
import sys

sys.path.append("/SNS/REF_L/shared/reduction")
from lr_reduction import time_resolved

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True)

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Time-resolved at 30Hz
    dynanic30_parser = subparsers.add_parser("dynamic30Hz", help="Reduce time-resolved 30Hz [-h for help]")
    dynanic30_parser.add_argument("meas_run_30Hz", type=int, help="Run number for the data to be processed")
    dynanic30_parser.add_argument(
        "ref_run_30Hz",
        type=str,
        help="Run number for the reference 30Hz data, measured at the same settings as the data to be processed",
    )
    dynanic30_parser.add_argument("ref_data_60Hz", type=str, help="Reference R(Q), measured at 60Hz")
    dynanic30_parser.add_argument("template_30Hz", type=str, help="File path for the 30Hz reduction template")
    dynanic30_parser.add_argument("time_interval", type=float, help="Time interval to use, in seconds")
    dynanic30_parser.add_argument("output_dir", type=str, help="Output directory")
    dynanic30_parser.add_argument(
        "--scan_index", type=int, dest="scan_index", help="Template scan index", required=False, default=1
    )
    dynanic30_parser.add_argument("--no-plot", dest="create_plot", action="store_false")
    dynanic30_parser.set_defaults(create_plot=True)
    dynanic30_parser.add_argument("--qsumming", dest="q_summing", action="store_true")
    dynanic30_parser.set_defaults(q_summing=False)

    # Standard template reduction for time-resolved data
    dynanic_parser = subparsers.add_parser("template", help="Reduce time-resolved [-h for help]")
    dynanic_parser.add_argument("run", type=int, help="Run number for the data to be processed")
    dynanic_parser.add_argument("template", type=str, help="File path for the reduction template")
    dynanic_parser.add_argument("time_interval", type=float, help="Time interval to use, in seconds")
    dynanic_parser.add_argument("output_dir", type=str, help="Output directory")
    dynanic_parser.add_argument(
        "--scan_index", type=int, dest="scan_index", help="Template scan index", required=False, default=1
    )
    dynanic_parser.add_argument(
        "--offset", type=float, dest="theta_offset", help="Theta offset", required=False, default=None
    )
    dynanic_parser.add_argument("--no-plot", dest="create_plot", action="store_false")
    dynanic_parser.set_defaults(create_plot=True)
    dynanic_parser.add_argument("--qsumming", dest="q_summing", action="store_true")
    dynanic_parser.set_defaults(q_summing=False)

    # Parse arguments
    args = parser.parse_args()

    if args.command == "dynamic30Hz":
        print("Time-resolved reduction with D2O reference: run %s" % args.meas_run_30Hz)
        reduced = time_resolved.reduce_30Hz_slices(
            args.meas_run_30Hz,
            args.ref_run_30Hz,
            args.ref_data_60Hz,
            args.template_30Hz,
            time_interval=args.time_interval,
            output_dir=args.output_dir,
            scan_index=args.scan_index,
            create_plot=args.create_plot,
            q_summing=args.q_summing,
        )
    else:
        print("Time-resolved reduction with template: run %s" % args.run)
        reduced = time_resolved.reduce_slices(
            args.run,
            args.template,
            time_interval=args.time_interval,
            output_dir=args.output_dir,
            scan_index=args.scan_index,
            theta_offset=args.theta_offset,
            create_plot=args.create_plot,
        )
