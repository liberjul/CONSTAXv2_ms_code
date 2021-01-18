import xlsxwriter
import pandas as pd

db_line = ",," + ",".join(["UNITE"]*12) + "," + ",".join(["SILVA"]*12) + "\n"
reg_line = ",," + ",".join(["Full"]*4) + "," + ",".join(["ITS1"]*4) + "," + ",".join(["ITS2"]*4) + "," + ",".join(["Full"]*4) + "," + ",".join(["V3-4"]*4) + "," + ",".join(["V4"]*4) + "\n"
met_line = "Partition Level,Classifier," + ",".join(["EPQ", "MC", "OC", "sensitivity"]*6) + "\n"

with open("../../tables/region_partition_summary.csv", "r") as ifile, open("../../tables/region_partition_summary_formatted.csv", "w") as ofile:
    ofile.write(db_line)
    ofile.write(reg_line)
    ofile.write(met_line)
    line = ifile.readline()
    line = ifile.readline()
    while line != "":
        ofile.write(",".join(line.split(",")[1:]))
        line = ifile.readline()

workbook = xlsxwriter.Workbook('../../tables/region_partition_summary_formatted.xlsx')
worksheet = workbook.add_worksheet()

row = 0
col = 0
merge_format = workbook.add_format({
    'align': 'center',
    'valign': 'vcenter'})
head_format = workbook.add_format({
    'bottom': 1,
    'align': 'center',
    'valign': 'vcenter'})
dark_format = workbook.add_format({
#     'bottom': 1,
    'align': 'center',
    'valign': 'vcenter',
    'fg_color': '#A9A9A9'})
light_format = workbook.add_format({
#     'bottom': 1,
    'align': 'center',
    'valign': 'vcenter',
    'fg_color': '#D3D3D3'})
with open("../../tables/region_partition_summary_formatted.csv", "r") as ifile:
    line = ifile.readline()
    while line != "":
        spl_line = line.strip().split(",")
        for item in spl_line:
            item = item.strip('"')
            if item == "NA":
                item = ""
            if row == 2 and col > 1:
                worksheet.write(row, col, item, head_format)
            elif row > 2 and (row % 2 == 0) and col > 1 and (col % 2 != 0):
                worksheet.write(row, col, item, dark_format)
            elif row > 2 and (row % 2 == 0) and col > 0:
                worksheet.write(row, col, item, light_format)
            elif row > 2 and (row % 2 != 0) and col > 1 and (col % 2 != 0):
                worksheet.write(row, col, item, light_format)
            else:
                worksheet.write(row, col, item, merge_format)
            col += 1
        row += 1
        col = 0
        line = ifile.readline()

worksheet.merge_range('C1:N1', 'UNITE', merge_format)
worksheet.merge_range('O1:Z1', 'SILVA', merge_format)
worksheet.merge_range('C2:F2', 'Full', merge_format)
worksheet.merge_range('G2:J2', 'ITS1', merge_format)
worksheet.merge_range('K2:N2', 'ITS2', merge_format)
worksheet.merge_range('O2:R2', 'Full', merge_format)
worksheet.merge_range('S2:V2', 'V3-4', merge_format)
worksheet.merge_range('W2:Z2', 'V4', merge_format)

worksheet.merge_range('A4:A14', 'Family', merge_format)
worksheet.merge_range('A15:A25', 'Genus', merge_format)
worksheet.set_column('A:A', 12)
worksheet.set_column('B:B', 17)

workbook.close()
