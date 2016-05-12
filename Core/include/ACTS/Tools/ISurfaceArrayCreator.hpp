///////////////////////////////////////////////////////////////////
// ISurfaceArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ISURFACEARRAYCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ISURFACEARRAYCREATOR_H 1

// Geometry module
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"
// STL
#include <vector>

namespace Acts {

  /** forward declarations & typedef */
  class Surface;
  typedef BinnedArray<const Surface*> SurfaceArray;
  
  /** @class ISurfaceArrayCreator
  
    Interface class ISurfaceArrayCreators, it inherits from IAlgTool. 
        
    @author Andreas.Salzburger@cern.ch
  */
  
  class ISurfaceArrayCreator {
    
    public:
      /**Virtual destructor*/
      virtual ~ISurfaceArrayCreator(){}

      /** SurfaceArrayCreator interface method - create an array in a cylinder, binned in phi, z */
      virtual SurfaceArray* surfaceArrayOnCylinder(const std::vector<const Surface*>& surfaces,
                                                   double R, double minPhi, double maxPhi, double halfZ, 
                                                   size_t binsPhi, size_t binsZ, 
                                                   std::shared_ptr<Transform3D> transform = nullptr) const = 0; 

      /** SurfaceArrayCreator interface method - create an array on a disc, binned in r, phi */
      virtual SurfaceArray* surfaceArrayOnDisc(const std::vector<const Surface*>& surfaces,
                                               double rMin, double rMax, double minPhi, double maxPhi,
                                               size_t binsR, size_t binsZ,
                                               const std::vector<double>& rBoundaries = {},
                                               std::shared_ptr<Transform3D> transform = nullptr) const = 0; 

      /** SurfaceArrayCreator interface method - create an array on a plane */
      virtual SurfaceArray* surfaceArrayOnPlane(const std::vector<const Surface*>& surfaces,
                                                double halflengthX, double halflengthY, 
                                                size_t binsX, size_t binsY,
                                                std::shared_ptr<Transform3D> transform = nullptr) const = 0; 
   
  };


} // end of namespace


#endif // ACTS_GEOMETRYINTERFACES_ISURFACEARRAYCREATOR_H