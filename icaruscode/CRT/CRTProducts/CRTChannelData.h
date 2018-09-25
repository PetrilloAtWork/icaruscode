#ifndef CRTChannelData_h_
#define CRTChannelData_h_

#include <stdint.h>
#include <vector>

namespace icarus {
namespace crt {

  class CRTChannelData {

    public:
      CRTChannelData();
      CRTChannelData(uint32_t chan, int time0, int time1, uint32_t q, std::vector<int> trackid);
      virtual ~CRTChannelData();
      uint32_t Channel() const;
      int T0() const;
      int T1() const;
      uint32_t ADC() const;
      std::vector<int> TrackID() const;
      void SetADC(uint32_t adc);

    private:
      uint32_t fChannel;
      int fT0;
      int fT1;
      uint32_t fAdc;
      std::vector<int> fTrackID;
  };

 }
}

#endif