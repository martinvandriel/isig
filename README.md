# ISIG: Instaseis Sismospher Input Generator

Generate Sismospher Input:
```
python -m isig -h
usage: python -m isig [-h] -m MODEL_FILE [-p PERIOD] [-d MAX_DEPTH]
                      [--min_dist MIN_DIST] [--max_dist MAX_DIST]
                      [-e ELEMENTS_PER_WAVELENGTH] [-n NPOL]

Generate a Instaseis input for Sismosphere.

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL_FILE, --model_file MODEL_FILE
                        path to 1D model in deck file format (or any other
                        format compatible with salvus mesher) (default: None)
  -p PERIOD, --period PERIOD
                        Shortest period to resolve. (default: 50.0)
  -d MAX_DEPTH, --max_depth MAX_DEPTH
                        Maximum source depth in km. (default: 100.0)
  --min_dist MIN_DIST   Minimum epicentral distance in degrees. (default: 0.0)
  --max_dist MAX_DIST   Maximum epicentral distance in degrees. (default:
                        180.0)
  -e ELEMENTS_PER_WAVELENGTH, --elements_per_wavelength ELEMENTS_PER_WAVELENGTH
                        Number of Elements per Wavelength. (default: 2.0)
  -n NPOL, --npol NPOL  Polynomial order used for interpolation. (default: 4)

```
