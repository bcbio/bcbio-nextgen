#!/usr/bin/env python
"""Retrieve a high level summary report of sequencing done in a month.

Usage:
    sequencing_report.py --month=<month> --year=<year> <YAML post process config>

month and year are both optional, in which case we'll default to this month
in the current year, which is the standard report of interest.

A month runs from the start of the 15th of the previous month to the end of the
14th in the current month, and tracks all projects which had their states set to
complete in this period.
"""
import sys
import csv
import calendar
from datetime import datetime
from optparse import OptionParser

import yaml

from bcbio.galaxy.api import GalaxyApiAccess

def main(config_file, month, year):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    galaxy_api = GalaxyApiAccess(config["galaxy_url"],
        config["galaxy_api_key"])
    smonth, syear = (month - 1, year) if month > 1 else (12, year - 1)
    start_date = datetime(syear, smonth, 15, 0, 0, 0)
    # last day calculation useful if definition of month is
    # from first to last day instead of 15th-15th
    #(_, last_day) = calendar.monthrange(year, month)
    end_date = datetime(year, month, 14, 23, 59, 59)
    out_file = "%s_%s" % (start_date.strftime("%b"),
            end_date.strftime("%b-%Y-sequencing.csv"))
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow([
            "Date", "Product", "Payment", "Researcher", "Lab", "Email",
            "Project", "Sample", "Description", "Genome", "Flowcell",
            "Lane", "Notes"])
        for s in galaxy_api.sqn_report(start_date.isoformat(),
                end_date.isoformat()):
            f_parts = s["sqn_run"]["run_folder"].split("_")
            flowcell = "_".join([f_parts[0], f_parts[-1]])
            writer.writerow([
                s["sqn_run"]["date"],
                s["sqn_type"],
                s["project"]["payment_(fund_number)"],
                s["project"]["researcher"],
                s["project"]["lab_association"],
                s["project"]["email"],
                s["project"]["project_name"],
                s["name"],
                s["description"],
                s["genome_build"],
                flowcell,
                s["sqn_run"]["lane"],
                s["sqn_run"]["results_notes"]])

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-m", "--month", dest="month")
    parser.add_option("-y", "--year", dest="year")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        print __doc__
        sys.exit()
    cur_year, cur_month = datetime.now().timetuple()[:2]
    if not options.month:
        options.month = cur_month
    if not options.year:
        options.year = cur_year
    main(args[0], int(options.month), int(options.year))
