#pragma once

#include "image_io.hpp"

// One trace. Width and length run between its ridges (the bright-edge peaks of the image's row and
// column profiles); the box encloses the whole trace: its outermost bright edge pixels plus 3 px. Where a
// neck connects the horizontal trace to the vertical one, its box ends where its own end cap would end.
struct Trace {
    int id = 0;  // 1: horizontal trace (left), 2: vertical trace (right)
    bool vertical = false;
    double along0 = 0, along1 = 0;      // ridges along the trace: x for the horizontal one, y for the vertical one
    double across0 = 0, across1 = 0;    // ridges across it
    double alongOut0 = 0, alongOut1 = 0, acrossOut0 = 0, acrossOut1 = 0;  // the enclosing box
    bool connected = false;             // horizontal trace only: a neck joins it to the vertical trace
    double neckWidth = 0;               // the neck's narrowest width, ridge to ridge, px
    double neckAt = 0, neckFrom = 0, neckTo = 0;  // that narrowest cut (the waist): x, and its two edges' y
    double neckX0 = 0, neckX1 = -1, neckY0 = 0, neckY1 = -1;  // the neck's box (between box 1 and box 2)
    std::vector<double> neckProfile;    // the neck's cross-section at the waist: intensity at y = neckFrom0 + i over its box's rows
    double neckFrom0 = 0, neckEdge0 = 0, neckEdge1 = 0;  // neckEdge*: the waist's two edges, y
    double contrast = 0;                // the trace's peak contrast over the background's (below 2: no trace)
    std::vector<double> crossSection;   // mean along the trace (following its slight tilt) at across = crossFrom + i, box edge to box edge
    double crossFrom = 0;
    int plotFirst = 0, plotLast = 0;    // crossSection indices to plot

    double length() const { return along1 - along0; }
    double width() const { return across1 - across0; }
    double x0() const { return vertical ? acrossOut0 : alongOut0; }
    double x1() const { return vertical ? acrossOut1 : alongOut1; }
    double y0() const { return vertical ? alongOut0 : acrossOut0; }
    double y1() const { return vertical ? alongOut1 : acrossOut1; }
};

// How an image shows its traces: thin bright edges around a dark inside (as imaged before processing), or
// thick bright rims around an inside about as dark as the background (after processing). Each is analysed
// its own way.
enum class ImageStyle { Edges, Rims };

struct Analysis {
    ImageStyle style = ImageStyle::Edges;
    std::vector<Trace> traces;  // the horizontal trace (if found), then the vertical one (if found)
    std::string error;          // a manufacture error, e.g. "manufacture error (dot array)"; empty if none
};

// Expects one horizontal trace entering from the left edge and one vertical trace to its right.
// A trace that is not there is left out.
Analysis analyzeImage(const Gray& img);
