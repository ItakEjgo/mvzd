import osmium
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Extract OSM Full History PBF into yearly CSVs.')
parser.add_argument('--input', required=True, help='Path to .osh.pbf file')
parser.add_argument('--outdir', required=True, help='Directory to store yearly CSVs')
args = parser.parse_args()

OUTPUT_DIR = args.outdir
os.makedirs(OUTPUT_DIR, exist_ok=True)

class YearlyExtractor(osmium.SimpleHandler):
    def __init__(self):
        super(YearlyExtractor, self).__init__()
        self.files = {}
        self.count = 0

    def get_file(self, year):
        if year not in self.files:
            path = os.path.join(OUTPUT_DIR, f"{year}.csv")
            f = open(path, 'w', encoding='utf-8')
            f.write("node_id,version,visible,changeset,timestamp,uid,lon,lat\n")
            self.files[year] = f
        return self.files[year]

    def node(self, n):
        self.count += 1
        
        year = n.timestamp.year
        f = self.get_file(year)

        visible = not n.deleted
        
        lon = n.location.lon if n.location.valid() else ""
        lat = n.location.lat if n.location.valid() else ""
        ts_str = n.timestamp.isoformat()

        row = f"{n.id},{n.version},{visible},{n.changeset},{ts_str},{n.uid},{lon},{lat}\n"
        f.write(row)
        
    def close(self):
        for f in self.files.values():
            f.close()

if __name__ == '__main__':
    print(f"Extracting {args.input} to {OUTPUT_DIR}...")
    handler = YearlyExtractor()
    handler.apply_file(args.input)
    handler.close()
    print(f"Extraction complete. Total nodes processed: {handler.count}")
