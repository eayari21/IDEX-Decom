import requests
import csv
import xml.etree.ElementTree as ET
from struct import unpack

# --- Load ELID to Command Name Mapping from XTCE ---
def load_elid_mapping(xtce_file):
    elid_map = {}
    ns = {'xtce': "http://www.omg.org/space/xtce"}
    tree = ET.parse(xtce_file)
    root = tree.getroot()
    for enum in root.findall(".//xtce:EnumeratedParameterType[@name='IDX_EVT.ELID_EVTPKT']//xtce:Enumeration", ns):
        val = int(enum.attrib['value'])
        label = enum.attrib['label']
        elid_map[val] = label
    return elid_map

# --- Parse CCSDS Packets ---
def parse_ccsds_packets(binary_data):
    offset = 0
    packets = []
    while offset + 6 <= len(binary_data):
        header = binary_data[offset:offset+6]
        apid = int.from_bytes(header[0:2], "big") & 0x07FF
        version = (int.from_bytes(header[0:2], "big") >> 13) & 0x07
        sequence = int.from_bytes(header[2:4], "big")
        length = int.from_bytes(header[4:6], "big") + 1
        total_size = length + 6
        if offset + total_size > len(binary_data):
            break
        packets.append({
            "raw": binary_data[offset:offset + total_size],
            "apid": apid,
            "version": version,
            "sequence": sequence,
            "length": length
        })
        offset += total_size
    return packets

# --- Decode EVTLOG Payload ---
def decode_evtlog(packet_bytes, elid_map):
    payload = packet_bytes[6:]  # Skip 6-byte CCSDS header
    shcoarse = unpack(">I", payload[0:4])[0]
    shfine = unpack(">H", payload[4:6])[0]
    elsec = unpack(">I", payload[6:10])[0]
    elssec = unpack(">H", payload[10:12])[0]
    elslice = payload[12]
    elid = payload[13]
    label = elid_map.get(elid, f"UNKNOWN_{elid}")
    return shcoarse, shfine, elsec, elssec, elslice, elid, label

# --- Auth + Download Binary ---
username = "eayari"
password = "xukco8-jujpax-fAkxoz"
url = "https://lasp.colorado.edu/ops/imap/poda/dap2/packets/SID2/IDX_EVTLOG.bin?time>=1742438400000000&time<1743292800000000"
# url = "https://lasp.colorado.edu/ops/imap/poda/dap2/packets/SID2/IDX_EVTLOG.bin?time>=1262304000000000"

response = requests.get(url, auth=(username, password))
with open("idx_evtlog_2024_2025.bin", "wb") as f:
    f.write(response.content)

# --- Load XTCE Mapping ---
xtce_file = "idex_housekeeping_packet_definition.xml"  # Replace with actual path
elid_map = load_elid_mapping(xtce_file)
print(f"✅ Loaded {len(elid_map)} ELID mappings from XTCE")

# --- Parse and Decode Packets ---
with open("idx_evtlog_2024_2025.bin", "rb") as f:
    data = f.read()

parsed = parse_ccsds_packets(data)
print(f"✅ Parsed {len(parsed)} packets")

# --- Write Full CSV ---
with open("evtlog_full_decoded.csv", "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow([
        "SEQ", "VERSION", "LENGTH",
        "SHCOARSE", "SHFINE",
        "ELSEC_EVTPKT", "ELSSEC_EVTPKT",
        "ELSLICE_EVTPKT", "ELID_EVTPKT", "COMMAND_NAME"
    ])
    count = 0
    for pkt in parsed:
        if pkt["apid"] != 4:
            continue
        fields = decode_evtlog(pkt["raw"], elid_map)
        writer.writerow([
            pkt["sequence"], pkt["version"], pkt["length"],
            *fields
        ])
        count += 1

print(f"✅ Wrote {count} EVTLOG entries to evtlog_full_decoded.csv")
