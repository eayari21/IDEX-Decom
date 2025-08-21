# IDEX-Decom
A set of modules for processing data from the Interstellar Dust Experiment (IDEX) aboard NASA's Interstellar Mapping and Acceleration Probe (IMAP) mission.

IDEX-Decom provides a pipeline for processing Interstellar Dust Experiment (IDEX) telemetry from raw Level-0 packets to Level-2C science products.

## Repository layout

```
IDEX-Decom/
├── L0->L1A/
│   └── idex-decom/        # packet parsing and Rice decompression
├── L1A->L1B/
│   └── idex-decom/
├── L1B->L2A/
│   └── idex-decom/
├── L2A->L2B/
│   └── idex-decom/
├── L2B->L2C/
│   └── idex-decom/
└── vectorize_submodules.sh
```

Each stage converts data to a higher level product and contains scripts, sample `Data/`, optional `HDF5/` outputs, plotting utilities, and a local README describing usage.

## Key components
- **Packet handling** – `idex_packet.py` and `science_tool.py` parse XTCE-defined packet binaries and decompress Rice-coded waveforms.
- **Rice compression** – `rice_decode.py` and `david_rice_decode.py` implement 10- and 12-bit Rice/Golomb decoding.
- **OASIS integration** – `read_from_ois.py` connects to the OASIS socket server to record raw telemetry streams.
- **Support files** – `idex_combined_science_definition.xml` defines packet structure, dependency lists are in `requirements.txt` and `idex.txt`, and `vectorize_submodules.sh` turns submodules into normal directories.

## Getting started
1. **Set up an environment**
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install -r L0->L1A/idex-decom/idex.txt
   ```
   Ensure the `output_from_ois/` directory exists for captured telemetry.

2. **Run the Level-0 to Level-1A example**
   ```bash
   python L0->L1A/idex-decom/science_tool.py Data/example_packet.bin
   ```
   Explore equivalent scripts in subsequent folders to generate higher level products.

## Suggested next steps
- Review `idex_combined_science_definition.xml` and LASP's `lasp_packets` documentation for the XTCE packet format.
- Examine advanced tools such as `SNR_Calculator.py`, `time2mass.py`, and mapping utilities in `imap_processing`.
- Learn external libraries like `cdflib` for CDF output and the Rice compression algorithm used throughout the pipeline.

These resources should help newcomers contribute to or extend the IDEX data processing pipeline.
