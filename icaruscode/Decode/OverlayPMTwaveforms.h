// Full documented header implementation
#ifndef OVERLAYPMTWAVEFORMS_H
#define OVERLAYPMTWAVEFORMS_H

// Include necessary headers
#include <vector>
#include <iostream>

// OverlayPMTwaveforms class definition
class OverlayPMTwaveforms {
public:
    OverlayPMTwaveforms();
    ~OverlayPMTwaveforms();

    void overlay();
    std::vector<double> getWaveforms();

private:
    std::vector<double> waveforms;
};

#endif // OVERLAYPMTWAVEFORMS_H