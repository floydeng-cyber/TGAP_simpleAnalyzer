# in-Python port

## Python syntax
```python
import sys
sys.path.insert(0, "path/to/etchAnalyzer-hitachi/src")
from hitachi import width, gap, d_width, d_gap
```

| Function | Returns |
|---|---|
| `width(structure, wafer, chip, row, col, kind)` | number
| `gap(structure, wafer, chip, row, col, kind)` | number
| `d_width(structure, wafer, chip, row, col)` | Δpost-pre
| `d_gap(structure, wafer, chip, row, col)` | Δpost-pre

# Python example

```python
>>> width("TGAP", 1, "0,0", 4, 4)
255.51 nm
>>> width("TGAP", 1, "0,0", 4, 4, kind="post")
533.47 nm
>>> gap("TGAP", 1, "0,0", 4, 4)
370.97 nm
>>> d_width("TGAP", 1, "0,0", 4, 4)
277.96 nm
>>> d_gap("TGAP", 1, "0,0", 4, 4)
-174.9 nm
```


# csv-output port
```bash
python3 src/hitachi.py [STRUCTURE] [WAFER] [--kind pre|post] [--chip X,Y] [--row R] [--col C] [--change] [--fields]
```

| Parameter | Values |
|---|---|
| `STRUCTURE` | TGAP DTGAP TTGAP UL ETE TETE FLG TCL |
| `WAFER` | 1 2 3 4 5 6 7 8 9 |
| `--kind` | pre post | pre-etch, post-etch
| `--chip` | -4,0 0,-5 0,0 0,5 4,0 | chip_number (x,y)
| `--row` | 1 2 3 4 5 6 | row_in_chip
| `--col` | 1 2 3 4 5 6 | column_in_chip
| `--change` | --change | to show Δpost-pre
| `--fields` | width gap x_um y_um |

## Example

```
$ python3 src/hitachi.py TGAP 1 --chip 0,0 --row 4 --col 4
structure,wafer,kind,chip,grid_col,grid_row,x_um,y_um,width,gap
TGAP,1,pre,"0,0",4,4,4408.0,2706.0,255.51,370.97
TGAP,1,post,"0,0",4,4,4414.0,2693.0,533.47,196.07

$ python3 src/hitachi.py TGAP 1 --chip 0,0 --row 4 --col 4 --change
structure,wafer,chip,grid_col,grid_row,pre_width,post_width,d_width,pre_gap,post_gap,d_gap
TGAP,1,"0,0",4,4,255.51,533.47,277.96,370.97,196.07,-174.9
```