# trace_boxes

Finds the horizontal and the vertical trace in TIFF images, draws an axis-aligned box around each,
and plots each trace's cross-section. Plain C++17 and CMake, with no other dependencies.

## Build and run

You need a C++17 compiler and CMake 3.25 or newer:
- **macOS:** `xcode-select --install`, then `brew install cmake`.
- **Windows:** Visual Studio 2022 with "Desktop development with C++" (includes CMake).
- **Linux:** `build-essential`, `cmake`, and `zenity` for the file dialog.

**VS Code:** open the folder *holding* `floAnalyzer` (its `.vscode/` points CMake at `floAnalyzer/`)
and accept the recommended extensions (C/C++ and CMake Tools).
If asked, pick the `release` preset. Then run with any of these:
- the **▷ in the status bar** at the bottom (CMake Tools: builds, then runs);
- **F5** (builds, then runs under the debugger);
- the **▷ at the editor's top right**, if you have the Code Runner extension: the project's
  settings point it at the command below.

A run button that compiles only the open file (`g++ main.cpp`) cannot build this program,
because it spans several source files. The project's settings turn those buttons off or redirect them.

**Visual Studio:** File → Open → Folder, pick `trace_boxes.exe` as the startup item, then ▶.

**Terminal (any OS), from inside `floAnalyzer/`:**

```bash
cmake --workflow --preset run                   # configure, build and run: opens the file dialog
build/release/trace_boxes a.tif some/folder     # or pass files and/or folders directly
build/release/trace_boxes -o out some/folder
```

The dialog accepts TIFF files and folders together on macOS. On Windows and Linux the dialog
picks files only; pass folders on the command line instead. Folders are searched recursively.

## Output

Three kinds of output:
- **`<name>_traces.svg`** (open it in any browser), for every image with a neck.
- **`analysisResult.xlsx`**, holding the measurements of every image without a neck, and the list of
  special cases.
- **`similarityResult_<type>.svg`**, plotting every image's HT and VT similarity factors.
- **`withinInRange.xlsx`**, the images within the standard deviation range.

Each SVG holds:
- **Left:** the image with each trace's box, dashed long axis and id. Where the traces are
  connected, a dashed red box encloses the neck and a red line marks its narrowest cut.
- **Right:** one subplot per trace. The x-axis is the position across the trace (px from its left
  boundary); the y-axis is the intensity averaged along the trace (0–255). Dotted lines and dots
  mark the ridges, which are the boundaries. All subplots share one y-axis.
- **For a connection, one more subplot:** the neck's cross-section at the waist, i.e. intensity at
  each row of the neck box, averaged over 7 columns around the waist. Its dotted lines and dots mark
  the waist's two edges; its title gives the neck's narrowest width.

`analysisResult.xlsx` sits directly in `analysisOutput/`. It is written fresh on every run.
- **One sheet per image folder,** named after the folder. A name longer than Excel's 31 characters
  keeps its start and end around a `~`.
- **Each sheet:** a header row, then one row per image.
- **Columns:** the image name in three parts, then the horizontal trace's (HT) width, length and
  depth, then the vertical trace's (VT).
  - **Prefix:** up to the first `_`, e.g. `E06_`.
  - **Image:** from there to the next `-`, e.g. `M0002`.
  - **Suffix:** the rest, e.g. `-02MS.tif`.
  - Prefix and suffix are shown only where they change from the row above.
- **Width and length** are as measured, ridge to ridge, in px.
- **Depth** is the mean height of the two ridges in the cross-section minus its lowest point
  between them, in intensity (0–255).
- **A last sheet, "Special cases",** lists every image with a neck or without both traces.
  - **Columns:** folder (shown where it changes), the image name in the same three parts, and a
    comment.
  - **Comment:** `neck width: N px`, `manufacture error (dot array)`, `only vertical trace`,
    `only horizontal trace` or `no trace`.

`similarityResult_<type>.svg` sits directly in `analysisOutput/`, next to the workbook. The type is the
selected folder's name up to its first `_` (e.g. `TGAP` for `TGAP_Part1`), with one plot per type.
- **Title:** "<type> Analysis Result".
- **Dots:** every image without a neck, coloured by folder. Hover a dot for its image and values.
- **x-axis, HT similarity factor:** the horizontal trace's width × length × depth.
- **y-axis, VT similarity factor:** the same for the vertical trace.
- **Normalized:** both axes show (value − mean) / standard deviation over all of the type's images,
  so the mean is 0 (solid line) and mean ± SD is ±1 (dashed lines, shaded band).
- **Table (also the legend):** each folder's normalized mean ± SD, and all images' raw mean ± SD,
  to convert back.

`withinInRange.xlsx` sits next to them. It lists the images within the standard deviation range:
within mean ± SD on both similarity factors, i.e. inside the plot's shaded band (mean and SD taken
per type). Its format is the same as `analysisResult.xlsx`'s folder sheets: one sheet per folder, the
same columns. It has no "Special cases" sheet, because images with a neck have no similarity factors.
Its cells are links (blue, underlined):
- **Prefix** (only in rows where it is shown): opens the raw data folder holding that image.
- **Image:** opens the raw image itself.

The links are absolute `file://` paths to where the images were when the program ran, so they break
if the raw data is moved or the workbook is opened on another computer.

Each run replaces the previous one's output: before writing, it deletes the earlier `*_traces.svg` and
`.csv` files in the output folders it writes to, and earlier `similarityResult*.svg` plots and
`withinInRange.xlsx` (nothing else).

Where outputs go, wherever the images are: `analysisOutput/` in the folder holding `floAnalyzer`.
The program finds that folder when it runs, from its own location (it walks up from the program to
`floAnalyzer`), so `floAnalyzer` can be moved anywhere and still works after a rebuild. The file dialog
starts in that folder too. If the program runs outside `floAnalyzer`, both use the current folder.
- An image selected directly: `analysisOutput/<image's folder name>/`.
- Images inside a selected folder: `analysisOutput/<folder name>/<same sub-folders>/`.
- With `-o DIR`, `DIR` replaces `analysisOutput`.

## Algorithm

### Two kinds of image

The same structures look different before and after processing, so each image is first sorted into
one of two kinds and then analysed its own way:
- **Edged** (as imaged before processing): thin bright edges around a dark trace, on a lighter
  background.
- **Rimmed** (after processing): thick, soft bright rims around a trace whose inside is about as
  dark as the background.

The sorting needs no metadata. It takes the pixels at least 4σ brighter than the background and
measures how long their runs along the image's rows are (the median run). Thin edges give runs of a
few px; thick rims give runs of tens of px. Images with a median run of 12 px or more are rimmed.
On the TGAP data this separates the two sets completely: January images reach at most 7 px, March
images at least 17 px.

Both analyses report the same measurements (boxes, width, length, depth, neck), so the outputs are
the same for either kind. Each SVG's title says which kind the image was taken for.

### Edged images

The program assumes the layout of these test images: one horizontal trace (id 1) entering from
the left edge, and one vertical trace (id 2) to its right. Boxes are axis-aligned and enclose the
whole trace. Width and length run between the trace's **ridges** (the peaks of its bright SEM edges).

Two maps are computed first. **Contrast** is how far each pixel (5×5 mean) is from the background,
in noise σ (background and σ come from the median and median deviation, re-estimated from pixels
within 3σ until stable). **Bright pixels** are 4σ above the background in a 3×3 mean, which keeps
thin edges and end caps.

For the vertical trace (columns of the whole image), then the horizontal trace (rows of the region
left of the vertical trace's box):

1. **Find it.** Average the contrast along the trace's direction into a profile across it. Its peak
   is the trace, and the stretch above half the peak is its rough band. If the peak is less than
   2× the background's level, the trace isn't there and is left out.
2. **Its ends.** Average the contrast over the band into a profile along the trace. From its peak,
   the dark body runs while the profile stays above half the body's **typical** level, i.e. the
   median of its strong slices (gaps up to 20 px). A single bright tip or speck would otherwise set
   the bar too high for a faint, patchy trace and cut it short. The run then continues while slices
   hold any bright pixel (gaps up to 10 px), which takes in the end caps. The horizontal trace
   always starts at the image's left edge, and its right end also stops at an **end cap** (a slice
   whose middle tenth is brighter than the background), so that the dark shading in the gap past a
   thick trace's tip isn't taken in.
3. **Neck.** The horizontal trace is connected to the vertical one by a neck if both hold:
   - **Close enough:** it ends within 1.5 trace widths + 20 px of the vertical trace.
   - **Continuous:** a continuous bright structure (the neck's edges) runs within its rows from
     inside its body to the vertical trace's edge. This is tested as a widest path on the **local
     contrast** (5×5 mean minus 21×21 mean, which cancels the dark haze around necks): of all
     paths across, the one whose dimmest pixel is brightest must stay above −1.35σ. Across a gap,
     the background between tip and vertical trace breaks every path.
   - **Its outline.** The background is flooded in from the image's border, through openings at
     least 15 px wide. It stops at:
     - **walls:** connected stretches of local contrast ≥ 1σ that hold at least 60 pixels of
       ≥ 2.5σ, so a faint stretch of a long wall still stops it while a noise speck doesn't;
     - **the inside of a trace:** a 13×13 mean more than 4σ below the background.

     So the background fills every bay up to the traces' outer walls but can't squeeze through a
     faint wall into a neck. In each column, the neck's outline runs from the trace's center line
     out to where the background starts. Bright texture inside the neck, and the trace's own cap
     curving in, are not the outline. Where the outline still jumps inward by more than 3 px from
     one column to the next (the background slipping past a faint wall), those columns are bridged
     by a straight line between the good columns on either side.
   - **Its width** runs ridge to ridge, like the traces' widths. A ridge lies as far inside the
     outline as the trace's own ridges lie inside its outline where it still has its full width.
     The **waist** is the narrowest column (median over 7), drawn in red.
   - **Thread necks.** Where the background passes right through a column, the neck is a single
     thread too thin for an outline. Its width there is where its local contrast is at least
     halfway from its surroundings to its peak.
   - **Box 1 ends where the trace's own end cap would end.** Going toward the neck, the outline's
     ridges first curve in, as for a free end cap, then turn outward into the neck. A free
     cap is a half-ellipse; its length is measured on the vertical trace's two free ends in the
     same image. Each column on the curving-in part (from where the trace drops below 90% of its
     width, to the inflection, i.e. the steepest narrowing) predicts the tip of such a cap through
     it. Box 1 ends at the median prediction plus 3 px. Its width and cross-section come from the
     full-width part.
   - **The neck box** runs from box 1 to box 2. Vertically it spans the rows that hold the neck's
     edges, including the arms flaring into the vertical trace: the outermost rows, within half a
     trace width + 20 px of the center, with at least 4 pixels of local contrast ≥ 3σ. It always
     includes the waist's two edges.
4. **Its sides.** The outermost row or column within 20 px of the rough band whose bright pixels
   cover at least 3% of the trace's length (and at least 2 px). That is the outside of the bright
   edge, even when the edge is slightly tilted.
5. **The box** is those extents plus 3 px.
6. **Ridges.** A fitted line through the trace's center lets the profiles follow a slight tilt, so
   the two edges aren't smeared together.
   - Across the trace (averaged along it): the brightest point on each side of the darkest point
     (found on a 5 px mean, then moved to the raw peak) gives the **width**.
   - Along its dark interior: the brightest point near each end (within one trace width) gives the
     **length**.
7. **Cross-section:** the across profile between the box edges.

Checked on all 893 TGAP_Part1 images:
- Every box encloses its trace, with no bright edge pixel left outside, apart from box 1 stopping at
  a neck.
- The 5 images missing a trace are on the "Special cases" sheet, along with the 75 necks.
- 75 connections are found. These are every E-series image at the designed neck positions
  (M0038/50/62, 110/122/134, 182/194/206, 254/266/278, 326/338/350) plus the 8 S-series images at
  those positions, and nothing else.
- 9 boxes have 9–15 px of slack at a faint, tapered tip.
- All 75 necks' outlines were checked by eye. In two (E08_M0278, E12_M0194) the top wall at the
  junction is almost invisible, and the background still creeps a little way into the neck there.

### Rimmed images

Filling in whatever a rim encloses fails where debris (with rims of its own) sticks to a trace. So
each trace is instead found from its own inside, walking outward:
1. **Rim pixels.** Pixels whose 5×5 mean is at least 4σ above the background *and* 40% of the way
   from the background to the rims' bright level (the image's 98th percentile). The second rule
   keeps out the dimmer film of debris some sites carry. Connected pieces shorter than 20% of the
   image are particles and are dropped.
2. **Manufacture errors.** If nothing but particles is left, and there are at least 3, the image is
   a dot array: listed as `manufacture error (dot array)`.
3. **Vertical trace.**
   - **Rims and centre line:** its two rims are the two fullest columns of rim pixels (each at least
     30% of the image). Its centre line runs between them. If no inside separates the rims (a thin
     trace), the middle of the fullest band.
   - **Sides:** from the centre line, each row is walked outward to the first rim and through it;
     the rim's outer edge is the outline there.
   - **Debris:** rows where the rim walked through is over 1.5× as thick as usual are left out,
     since there it runs into debris or the horizontal trace. The box side is the 10th percentile of
     the rest, from outside.
   - **Ends:** walked up and down the centre line, through the cap rims.
4. **Horizontal trace.**
   - **Columns:** left of the vertical trace, its two rims are the two fullest rows of rim pixels.
     From its centre line, each column is walked up and down the same way, from the left edge on,
     until no rim is found within two trace widths (past its end cap).
   - **Connected:** if its outline runs on into the vertical trace, it is connected. Its box then
     ends at the narrowest column within two widths of the vertical trace (the waist), and the rest
     is the neck.
5. **Width, length, depth and neck width.** Measured ridge to ridge on intensity profiles, as for
   edged images; the ridges are the rims' peaks.

Checked on the March W2 set (606 images):
- **Sorted:** 605 images were taken as rimmed; the 13 dot arrays among them are listed as
  manufacture errors. The program also reads the January set exactly as before, with none of its
  893 images taken for rimmed.
- **Necks:** 46 connections, all at E-series sites. Their boxes hug the traces' own rims, not the
  debris.
- **Missing traces:** 7 rimmed images have only a vertical trace (a stub or nothing on the left).
  One more, E06_M0134, is a blank field with two printed marks; it is taken as edged and also comes
  out as only vertical.

Supported TIFF files: first page, 8/16-bit, gray or RGB(A), uncompressed or LZW.
