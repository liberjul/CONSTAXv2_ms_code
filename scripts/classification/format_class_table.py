import xlsxwriter
import pandas as pd

db_line = ",," + ",".join(["UNITE"]*9) + "," + ",".join(["SILVA"]*9) + "\n"
reg_line = ",," + ",".join(["Full"]*3) + "," + ",".join(["ITS1"]*3) + "," + ",".join(["ITS2"]*3) + "," + ",".join(["Full"]*3) + "," + ",".join(["V3-4"]*3) + "," + ",".join(["V4"]*3) + "\n"
met_line = "Partition Level,Classifier," + ",".join(["EPQ", "MC", "OC"]*6) + "\n"

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

worksheet.merge_range('C1:K1', 'UNITE', merge_format)
worksheet.merge_range('L1:T1', 'SILVA', merge_format)
worksheet.merge_range('C2:E2', 'Full', merge_format)
worksheet.merge_range('F2:H2', 'ITS1', merge_format)
worksheet.merge_range('I2:K2', 'ITS2', merge_format)
worksheet.merge_range('L2:N2', 'Full', merge_format)
worksheet.merge_range('O2:Q2', 'V3-4', merge_format)
worksheet.merge_range('R2:T2', 'V4', merge_format)

worksheet.merge_range('A4:A16', 'Family', merge_format)
worksheet.merge_range('A17:A29', 'Genus', merge_format)
worksheet.set_column('A:A', 12)
worksheet.set_column('B:B', 17)
worksheet.set_column('C:T', 12)

workbook.close()
