import csv

class BedParser(object):

    def __init__(self, filename):
        self.filename = filename

    def intervals_generator(self):
        with open(self.filename, 'r') as fh:
            reader = csv.reader(fh, delimiter='\t')
            for row in reader:
                if len(row) >= 3 and row[0] not in ['track', 'browser']:
                    contig = str(row[0])
                    start = int(row[1])
                    end = int(row[2])
                    yield (contig, start, end)
