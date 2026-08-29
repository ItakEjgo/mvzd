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

class YearlyExtractor:
    def __init__(self, outdir):
        self.outdir = outdir
        self.files = {}
        self.buffers = {}
        self.count = 0

    def get_file_and_buffer(self, year):
        if year not in self.files:
            path = os.path.join(self.outdir, f"{year}.csv")
            f = open(path, 'w', encoding='utf-8')
            f.write("node_id,version,visible,changeset,timestamp,uid,lon,lat\n")
            self.files[year] = f
            self.buffers[year] = []
        return self.files[year], self.buffers[year]

    def process_file(self, input_file):
        # 尝试使用 FileProcessor（新版 osmium 提供），支持 Python 生成器，能完美响应 Ctrl+C
        if hasattr(osmium, 'FileProcessor'):
            self._process_with_fileprocessor(input_file)
        else:
            self._process_with_handler(input_file)

    def _process_with_fileprocessor(self, input_file):
        fp = osmium.FileProcessor(input_file, osmium.osm.osm_entity_bits.NODE)
        for n in fp:
            self._process_node(n)

    def _process_with_handler(self, input_file):
        print("⚠️ Warning: Your 'osmium' version is too old (missing FileProcessor).")
        print("⚠️ Please run: pip install --upgrade osmium (for much faster & safer parsing).")
        print("⚠️ Falling back to SimpleHandler. Ctrl+C might be unresponsive.")
        
        class FallbackHandler(osmium.SimpleHandler):
            def __init__(self, extractor):
                super().__init__()
                self.extractor = extractor
            def node(self, n):
                self.extractor._process_node(n)
            # Define empty handlers to force C++ to yield to Python, avoiding silent C++ blocking
            def way(self, w): pass
            def relation(self, r): pass
            
        handler = FallbackHandler(self)
        handler.apply_file(input_file)

    def _process_node(self, n):
        self.count += 1
        
        year = n.timestamp.year
        f, buf = self.get_file_and_buffer(year)

        visible = not n.deleted
        
        if self.count % 1000000 == 0:
            print(f"Processed {self.count} nodes...", flush=True)

        # Shift coordinates to strictly positive domain (first quadrant)
        lon = n.location.lon + 180.0 if n.location.valid() else ""
        lat = n.location.lat + 90.0 if n.location.valid() else ""
        ts_str = n.timestamp.isoformat()
        row = f"{n.id},{n.version},{visible},{n.changeset},{ts_str},{n.uid},{lon},{lat}\n"
        buf.append(row)
        if len(buf) >= 20000:
            f.write("".join(buf))
            buf.clear()
    def close(self):
        for year, f in self.files.items():
            buf = self.buffers[year]
            if buf:
                f.write("".join(buf))
                buf.clear()
            f.close()

if __name__ == '__main__':
    print(f"Extracting {args.input} to {OUTPUT_DIR}...")
    extractor = YearlyExtractor(OUTPUT_DIR)
    extractor.process_file(args.input)
    extractor.close()
    print(f"Extraction complete. Total nodes processed: {extractor.count}")
