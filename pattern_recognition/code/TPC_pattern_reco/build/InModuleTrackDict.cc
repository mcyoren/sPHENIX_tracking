// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME InModuleTrackDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "../InModuleTrack.h"
#include "../InModuleTrackv1.h"
#include "../InModuleTrackContainer.h"
#include "../InModuleTrackContainerv1.h"
#include "../FullTrack.h"
#include "../FullTrackv1.h"
#include "../FullTrackContainer.h"
#include "../FullTrackContainerv1.h"
#include "../FullTrackVertex.h"
#include "../FullTrackVertexv1.h"
#include "../FullTrackVertexContainer.h"
#include "../FullTrackVertexContainerv1.h"
#include "../FinalTrack.h"
#include "../FinalTrackv1.h"
#include "../FinalTrackContainer.h"
#include "../FinalTrackContainerv1.h"
#include "../FinalTrackVertex.h"
#include "../FinalTrackVertexv1.h"
#include "../FinalTrackVertexContainer.h"
#include "../FinalTrackVertexContainerv1.h"
#include "../TpcPolyTrack.h"
#include "../TpcPolyTrackv1.h"
#include "../TpcPolyTrackContainer.h"
#include "../TpcPolyTrackContainerv1.h"
#include "../TpcPolyCluster.h"
#include "../TpcPolyClusterv1.h"
#include "../TpcPolyClusterContainer.h"
#include "../TpcPolyClusterContainerv1.h"
#include "../TpcPolyClusterTrack.h"
#include "../TpcPolyClusterTrackv1.h"
#include "../TpcPolyClusterTrackContainer.h"
#include "../TpcPolyClusterTrackContainerv1.h"
#include "../Fitter.h"
#include "../IdealPadMap.h"
#include "../TpcPadMap.h"
#include "../TpcPadMapv1.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace Fitter {
   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *Fitter_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("Fitter", 0 /*version*/, "Fitter.h", 12,
                     ::ROOT::Internal::DefineBehavior((void*)nullptr,(void*)nullptr),
                     &Fitter_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *Fitter_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static void *new_InModuleTrack(void *p = nullptr);
   static void *newArray_InModuleTrack(Long_t size, void *p);
   static void delete_InModuleTrack(void *p);
   static void deleteArray_InModuleTrack(void *p);
   static void destruct_InModuleTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InModuleTrack*)
   {
      ::InModuleTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::InModuleTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("InModuleTrack", ::InModuleTrack::Class_Version(), "InModuleTrack.h", 15,
                  typeid(::InModuleTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::InModuleTrack::Dictionary, isa_proxy, 4,
                  sizeof(::InModuleTrack) );
      instance.SetNew(&new_InModuleTrack);
      instance.SetNewArray(&newArray_InModuleTrack);
      instance.SetDelete(&delete_InModuleTrack);
      instance.SetDeleteArray(&deleteArray_InModuleTrack);
      instance.SetDestructor(&destruct_InModuleTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InModuleTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::InModuleTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::InModuleTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_InModuleTrackv1(void *p = nullptr);
   static void *newArray_InModuleTrackv1(Long_t size, void *p);
   static void delete_InModuleTrackv1(void *p);
   static void deleteArray_InModuleTrackv1(void *p);
   static void destruct_InModuleTrackv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InModuleTrackv1*)
   {
      ::InModuleTrackv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::InModuleTrackv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("InModuleTrackv1", ::InModuleTrackv1::Class_Version(), "InModuleTrackv1.h", 9,
                  typeid(::InModuleTrackv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::InModuleTrackv1::Dictionary, isa_proxy, 4,
                  sizeof(::InModuleTrackv1) );
      instance.SetNew(&new_InModuleTrackv1);
      instance.SetNewArray(&newArray_InModuleTrackv1);
      instance.SetDelete(&delete_InModuleTrackv1);
      instance.SetDeleteArray(&deleteArray_InModuleTrackv1);
      instance.SetDestructor(&destruct_InModuleTrackv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InModuleTrackv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::InModuleTrackv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::InModuleTrackv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_InModuleTrackContainer(void *p = nullptr);
   static void *newArray_InModuleTrackContainer(Long_t size, void *p);
   static void delete_InModuleTrackContainer(void *p);
   static void deleteArray_InModuleTrackContainer(void *p);
   static void destruct_InModuleTrackContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InModuleTrackContainer*)
   {
      ::InModuleTrackContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::InModuleTrackContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("InModuleTrackContainer", ::InModuleTrackContainer::Class_Version(), "InModuleTrackContainer.h", 11,
                  typeid(::InModuleTrackContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::InModuleTrackContainer::Dictionary, isa_proxy, 4,
                  sizeof(::InModuleTrackContainer) );
      instance.SetNew(&new_InModuleTrackContainer);
      instance.SetNewArray(&newArray_InModuleTrackContainer);
      instance.SetDelete(&delete_InModuleTrackContainer);
      instance.SetDeleteArray(&deleteArray_InModuleTrackContainer);
      instance.SetDestructor(&destruct_InModuleTrackContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InModuleTrackContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::InModuleTrackContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::InModuleTrackContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_InModuleTrackContainerv1(void *p = nullptr);
   static void *newArray_InModuleTrackContainerv1(Long_t size, void *p);
   static void delete_InModuleTrackContainerv1(void *p);
   static void deleteArray_InModuleTrackContainerv1(void *p);
   static void destruct_InModuleTrackContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InModuleTrackContainerv1*)
   {
      ::InModuleTrackContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::InModuleTrackContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("InModuleTrackContainerv1", ::InModuleTrackContainerv1::Class_Version(), "InModuleTrackContainerv1.h", 10,
                  typeid(::InModuleTrackContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::InModuleTrackContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::InModuleTrackContainerv1) );
      instance.SetNew(&new_InModuleTrackContainerv1);
      instance.SetNewArray(&newArray_InModuleTrackContainerv1);
      instance.SetDelete(&delete_InModuleTrackContainerv1);
      instance.SetDeleteArray(&deleteArray_InModuleTrackContainerv1);
      instance.SetDestructor(&destruct_InModuleTrackContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InModuleTrackContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::InModuleTrackContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::InModuleTrackContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrack(void *p = nullptr);
   static void *newArray_FullTrack(Long_t size, void *p);
   static void delete_FullTrack(void *p);
   static void deleteArray_FullTrack(void *p);
   static void destruct_FullTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrack*)
   {
      ::FullTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrack", ::FullTrack::Class_Version(), "FullTrack.h", 14,
                  typeid(::FullTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrack::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrack) );
      instance.SetNew(&new_FullTrack);
      instance.SetNewArray(&newArray_FullTrack);
      instance.SetDelete(&delete_FullTrack);
      instance.SetDeleteArray(&deleteArray_FullTrack);
      instance.SetDestructor(&destruct_FullTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackv1(void *p = nullptr);
   static void *newArray_FullTrackv1(Long_t size, void *p);
   static void delete_FullTrackv1(void *p);
   static void deleteArray_FullTrackv1(void *p);
   static void destruct_FullTrackv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackv1*)
   {
      ::FullTrackv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackv1", ::FullTrackv1::Class_Version(), "FullTrackv1.h", 8,
                  typeid(::FullTrackv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackv1::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackv1) );
      instance.SetNew(&new_FullTrackv1);
      instance.SetNewArray(&newArray_FullTrackv1);
      instance.SetDelete(&delete_FullTrackv1);
      instance.SetDeleteArray(&deleteArray_FullTrackv1);
      instance.SetDestructor(&destruct_FullTrackv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackContainer(void *p = nullptr);
   static void *newArray_FullTrackContainer(Long_t size, void *p);
   static void delete_FullTrackContainer(void *p);
   static void deleteArray_FullTrackContainer(void *p);
   static void destruct_FullTrackContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackContainer*)
   {
      ::FullTrackContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackContainer", ::FullTrackContainer::Class_Version(), "FullTrackContainer.h", 9,
                  typeid(::FullTrackContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackContainer::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackContainer) );
      instance.SetNew(&new_FullTrackContainer);
      instance.SetNewArray(&newArray_FullTrackContainer);
      instance.SetDelete(&delete_FullTrackContainer);
      instance.SetDeleteArray(&deleteArray_FullTrackContainer);
      instance.SetDestructor(&destruct_FullTrackContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackContainerv1(void *p = nullptr);
   static void *newArray_FullTrackContainerv1(Long_t size, void *p);
   static void delete_FullTrackContainerv1(void *p);
   static void deleteArray_FullTrackContainerv1(void *p);
   static void destruct_FullTrackContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackContainerv1*)
   {
      ::FullTrackContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackContainerv1", ::FullTrackContainerv1::Class_Version(), "FullTrackContainerv1.h", 10,
                  typeid(::FullTrackContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackContainerv1) );
      instance.SetNew(&new_FullTrackContainerv1);
      instance.SetNewArray(&newArray_FullTrackContainerv1);
      instance.SetDelete(&delete_FullTrackContainerv1);
      instance.SetDeleteArray(&deleteArray_FullTrackContainerv1);
      instance.SetDestructor(&destruct_FullTrackContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackVertex(void *p = nullptr);
   static void *newArray_FullTrackVertex(Long_t size, void *p);
   static void delete_FullTrackVertex(void *p);
   static void deleteArray_FullTrackVertex(void *p);
   static void destruct_FullTrackVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackVertex*)
   {
      ::FullTrackVertex *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackVertex >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackVertex", ::FullTrackVertex::Class_Version(), "FullTrackVertex.h", 7,
                  typeid(::FullTrackVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackVertex::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackVertex) );
      instance.SetNew(&new_FullTrackVertex);
      instance.SetNewArray(&newArray_FullTrackVertex);
      instance.SetDelete(&delete_FullTrackVertex);
      instance.SetDeleteArray(&deleteArray_FullTrackVertex);
      instance.SetDestructor(&destruct_FullTrackVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackVertex*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackVertex*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackVertex*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackVertexv1(void *p = nullptr);
   static void *newArray_FullTrackVertexv1(Long_t size, void *p);
   static void delete_FullTrackVertexv1(void *p);
   static void deleteArray_FullTrackVertexv1(void *p);
   static void destruct_FullTrackVertexv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackVertexv1*)
   {
      ::FullTrackVertexv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackVertexv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackVertexv1", ::FullTrackVertexv1::Class_Version(), "FullTrackVertexv1.h", 7,
                  typeid(::FullTrackVertexv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackVertexv1::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackVertexv1) );
      instance.SetNew(&new_FullTrackVertexv1);
      instance.SetNewArray(&newArray_FullTrackVertexv1);
      instance.SetDelete(&delete_FullTrackVertexv1);
      instance.SetDeleteArray(&deleteArray_FullTrackVertexv1);
      instance.SetDestructor(&destruct_FullTrackVertexv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackVertexv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackVertexv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackVertexv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackVertexContainer(void *p = nullptr);
   static void *newArray_FullTrackVertexContainer(Long_t size, void *p);
   static void delete_FullTrackVertexContainer(void *p);
   static void deleteArray_FullTrackVertexContainer(void *p);
   static void destruct_FullTrackVertexContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackVertexContainer*)
   {
      ::FullTrackVertexContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackVertexContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackVertexContainer", ::FullTrackVertexContainer::Class_Version(), "FullTrackVertexContainer.h", 9,
                  typeid(::FullTrackVertexContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackVertexContainer::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackVertexContainer) );
      instance.SetNew(&new_FullTrackVertexContainer);
      instance.SetNewArray(&newArray_FullTrackVertexContainer);
      instance.SetDelete(&delete_FullTrackVertexContainer);
      instance.SetDeleteArray(&deleteArray_FullTrackVertexContainer);
      instance.SetDestructor(&destruct_FullTrackVertexContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackVertexContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackVertexContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackVertexContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FullTrackVertexContainerv1(void *p = nullptr);
   static void *newArray_FullTrackVertexContainerv1(Long_t size, void *p);
   static void delete_FullTrackVertexContainerv1(void *p);
   static void deleteArray_FullTrackVertexContainerv1(void *p);
   static void destruct_FullTrackVertexContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FullTrackVertexContainerv1*)
   {
      ::FullTrackVertexContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FullTrackVertexContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FullTrackVertexContainerv1", ::FullTrackVertexContainerv1::Class_Version(), "FullTrackVertexContainerv1.h", 10,
                  typeid(::FullTrackVertexContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FullTrackVertexContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::FullTrackVertexContainerv1) );
      instance.SetNew(&new_FullTrackVertexContainerv1);
      instance.SetNewArray(&newArray_FullTrackVertexContainerv1);
      instance.SetDelete(&delete_FullTrackVertexContainerv1);
      instance.SetDeleteArray(&deleteArray_FullTrackVertexContainerv1);
      instance.SetDestructor(&destruct_FullTrackVertexContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FullTrackVertexContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FullTrackVertexContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FullTrackVertexContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrack(void *p = nullptr);
   static void *newArray_FinalTrack(Long_t size, void *p);
   static void delete_FinalTrack(void *p);
   static void deleteArray_FinalTrack(void *p);
   static void destruct_FinalTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrack*)
   {
      ::FinalTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrack", ::FinalTrack::Class_Version(), "FinalTrack.h", 7,
                  typeid(::FinalTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrack::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrack) );
      instance.SetNew(&new_FinalTrack);
      instance.SetNewArray(&newArray_FinalTrack);
      instance.SetDelete(&delete_FinalTrack);
      instance.SetDeleteArray(&deleteArray_FinalTrack);
      instance.SetDestructor(&destruct_FinalTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackv1(void *p = nullptr);
   static void *newArray_FinalTrackv1(Long_t size, void *p);
   static void delete_FinalTrackv1(void *p);
   static void deleteArray_FinalTrackv1(void *p);
   static void destruct_FinalTrackv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackv1*)
   {
      ::FinalTrackv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackv1", ::FinalTrackv1::Class_Version(), "FinalTrackv1.h", 8,
                  typeid(::FinalTrackv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackv1::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackv1) );
      instance.SetNew(&new_FinalTrackv1);
      instance.SetNewArray(&newArray_FinalTrackv1);
      instance.SetDelete(&delete_FinalTrackv1);
      instance.SetDeleteArray(&deleteArray_FinalTrackv1);
      instance.SetDestructor(&destruct_FinalTrackv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackContainer(void *p = nullptr);
   static void *newArray_FinalTrackContainer(Long_t size, void *p);
   static void delete_FinalTrackContainer(void *p);
   static void deleteArray_FinalTrackContainer(void *p);
   static void destruct_FinalTrackContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackContainer*)
   {
      ::FinalTrackContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackContainer", ::FinalTrackContainer::Class_Version(), "FinalTrackContainer.h", 9,
                  typeid(::FinalTrackContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackContainer::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackContainer) );
      instance.SetNew(&new_FinalTrackContainer);
      instance.SetNewArray(&newArray_FinalTrackContainer);
      instance.SetDelete(&delete_FinalTrackContainer);
      instance.SetDeleteArray(&deleteArray_FinalTrackContainer);
      instance.SetDestructor(&destruct_FinalTrackContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackContainerv1(void *p = nullptr);
   static void *newArray_FinalTrackContainerv1(Long_t size, void *p);
   static void delete_FinalTrackContainerv1(void *p);
   static void deleteArray_FinalTrackContainerv1(void *p);
   static void destruct_FinalTrackContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackContainerv1*)
   {
      ::FinalTrackContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackContainerv1", ::FinalTrackContainerv1::Class_Version(), "FinalTrackContainerv1.h", 10,
                  typeid(::FinalTrackContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackContainerv1) );
      instance.SetNew(&new_FinalTrackContainerv1);
      instance.SetNewArray(&newArray_FinalTrackContainerv1);
      instance.SetDelete(&delete_FinalTrackContainerv1);
      instance.SetDeleteArray(&deleteArray_FinalTrackContainerv1);
      instance.SetDestructor(&destruct_FinalTrackContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackVertex(void *p = nullptr);
   static void *newArray_FinalTrackVertex(Long_t size, void *p);
   static void delete_FinalTrackVertex(void *p);
   static void deleteArray_FinalTrackVertex(void *p);
   static void destruct_FinalTrackVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackVertex*)
   {
      ::FinalTrackVertex *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackVertex >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackVertex", ::FinalTrackVertex::Class_Version(), "FinalTrackVertex.h", 7,
                  typeid(::FinalTrackVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackVertex::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackVertex) );
      instance.SetNew(&new_FinalTrackVertex);
      instance.SetNewArray(&newArray_FinalTrackVertex);
      instance.SetDelete(&delete_FinalTrackVertex);
      instance.SetDeleteArray(&deleteArray_FinalTrackVertex);
      instance.SetDestructor(&destruct_FinalTrackVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackVertex*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackVertex*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackVertex*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackVertexv1(void *p = nullptr);
   static void *newArray_FinalTrackVertexv1(Long_t size, void *p);
   static void delete_FinalTrackVertexv1(void *p);
   static void deleteArray_FinalTrackVertexv1(void *p);
   static void destruct_FinalTrackVertexv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackVertexv1*)
   {
      ::FinalTrackVertexv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackVertexv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackVertexv1", ::FinalTrackVertexv1::Class_Version(), "FinalTrackVertexv1.h", 7,
                  typeid(::FinalTrackVertexv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackVertexv1::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackVertexv1) );
      instance.SetNew(&new_FinalTrackVertexv1);
      instance.SetNewArray(&newArray_FinalTrackVertexv1);
      instance.SetDelete(&delete_FinalTrackVertexv1);
      instance.SetDeleteArray(&deleteArray_FinalTrackVertexv1);
      instance.SetDestructor(&destruct_FinalTrackVertexv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackVertexv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackVertexv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackVertexv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackVertexContainer(void *p = nullptr);
   static void *newArray_FinalTrackVertexContainer(Long_t size, void *p);
   static void delete_FinalTrackVertexContainer(void *p);
   static void deleteArray_FinalTrackVertexContainer(void *p);
   static void destruct_FinalTrackVertexContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackVertexContainer*)
   {
      ::FinalTrackVertexContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackVertexContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackVertexContainer", ::FinalTrackVertexContainer::Class_Version(), "FinalTrackVertexContainer.h", 9,
                  typeid(::FinalTrackVertexContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackVertexContainer::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackVertexContainer) );
      instance.SetNew(&new_FinalTrackVertexContainer);
      instance.SetNewArray(&newArray_FinalTrackVertexContainer);
      instance.SetDelete(&delete_FinalTrackVertexContainer);
      instance.SetDeleteArray(&deleteArray_FinalTrackVertexContainer);
      instance.SetDestructor(&destruct_FinalTrackVertexContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackVertexContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackVertexContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackVertexContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FinalTrackVertexContainerv1(void *p = nullptr);
   static void *newArray_FinalTrackVertexContainerv1(Long_t size, void *p);
   static void delete_FinalTrackVertexContainerv1(void *p);
   static void deleteArray_FinalTrackVertexContainerv1(void *p);
   static void destruct_FinalTrackVertexContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FinalTrackVertexContainerv1*)
   {
      ::FinalTrackVertexContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FinalTrackVertexContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("FinalTrackVertexContainerv1", ::FinalTrackVertexContainerv1::Class_Version(), "FinalTrackVertexContainerv1.h", 10,
                  typeid(::FinalTrackVertexContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FinalTrackVertexContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::FinalTrackVertexContainerv1) );
      instance.SetNew(&new_FinalTrackVertexContainerv1);
      instance.SetNewArray(&newArray_FinalTrackVertexContainerv1);
      instance.SetDelete(&delete_FinalTrackVertexContainerv1);
      instance.SetDeleteArray(&deleteArray_FinalTrackVertexContainerv1);
      instance.SetDestructor(&destruct_FinalTrackVertexContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FinalTrackVertexContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::FinalTrackVertexContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FinalTrackVertexContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyTrack(void *p = nullptr);
   static void *newArray_TpcPolyTrack(Long_t size, void *p);
   static void delete_TpcPolyTrack(void *p);
   static void deleteArray_TpcPolyTrack(void *p);
   static void destruct_TpcPolyTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyTrack*)
   {
      ::TpcPolyTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyTrack", ::TpcPolyTrack::Class_Version(), "TpcPolyTrack.h", 10,
                  typeid(::TpcPolyTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyTrack) );
      instance.SetNew(&new_TpcPolyTrack);
      instance.SetNewArray(&newArray_TpcPolyTrack);
      instance.SetDelete(&delete_TpcPolyTrack);
      instance.SetDeleteArray(&deleteArray_TpcPolyTrack);
      instance.SetDestructor(&destruct_TpcPolyTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyTrackv1(void *p = nullptr);
   static void *newArray_TpcPolyTrackv1(Long_t size, void *p);
   static void delete_TpcPolyTrackv1(void *p);
   static void deleteArray_TpcPolyTrackv1(void *p);
   static void destruct_TpcPolyTrackv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyTrackv1*)
   {
      ::TpcPolyTrackv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyTrackv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyTrackv1", ::TpcPolyTrackv1::Class_Version(), "TpcPolyTrackv1.h", 8,
                  typeid(::TpcPolyTrackv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyTrackv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyTrackv1) );
      instance.SetNew(&new_TpcPolyTrackv1);
      instance.SetNewArray(&newArray_TpcPolyTrackv1);
      instance.SetDelete(&delete_TpcPolyTrackv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyTrackv1);
      instance.SetDestructor(&destruct_TpcPolyTrackv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyTrackv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyTrackv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyTrackv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyTrackContainer(void *p = nullptr);
   static void *newArray_TpcPolyTrackContainer(Long_t size, void *p);
   static void delete_TpcPolyTrackContainer(void *p);
   static void deleteArray_TpcPolyTrackContainer(void *p);
   static void destruct_TpcPolyTrackContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyTrackContainer*)
   {
      ::TpcPolyTrackContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyTrackContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyTrackContainer", ::TpcPolyTrackContainer::Class_Version(), "TpcPolyTrackContainer.h", 9,
                  typeid(::TpcPolyTrackContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyTrackContainer::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyTrackContainer) );
      instance.SetNew(&new_TpcPolyTrackContainer);
      instance.SetNewArray(&newArray_TpcPolyTrackContainer);
      instance.SetDelete(&delete_TpcPolyTrackContainer);
      instance.SetDeleteArray(&deleteArray_TpcPolyTrackContainer);
      instance.SetDestructor(&destruct_TpcPolyTrackContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyTrackContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyTrackContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyTrackContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyTrackContainerv1(void *p = nullptr);
   static void *newArray_TpcPolyTrackContainerv1(Long_t size, void *p);
   static void delete_TpcPolyTrackContainerv1(void *p);
   static void deleteArray_TpcPolyTrackContainerv1(void *p);
   static void destruct_TpcPolyTrackContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyTrackContainerv1*)
   {
      ::TpcPolyTrackContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyTrackContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyTrackContainerv1", ::TpcPolyTrackContainerv1::Class_Version(), "TpcPolyTrackContainerv1.h", 10,
                  typeid(::TpcPolyTrackContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyTrackContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyTrackContainerv1) );
      instance.SetNew(&new_TpcPolyTrackContainerv1);
      instance.SetNewArray(&newArray_TpcPolyTrackContainerv1);
      instance.SetDelete(&delete_TpcPolyTrackContainerv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyTrackContainerv1);
      instance.SetDestructor(&destruct_TpcPolyTrackContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyTrackContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyTrackContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyTrackContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyCluster(void *p = nullptr);
   static void *newArray_TpcPolyCluster(Long_t size, void *p);
   static void delete_TpcPolyCluster(void *p);
   static void deleteArray_TpcPolyCluster(void *p);
   static void destruct_TpcPolyCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyCluster*)
   {
      ::TpcPolyCluster *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyCluster >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyCluster", ::TpcPolyCluster::Class_Version(), "TpcPolyCluster.h", 10,
                  typeid(::TpcPolyCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyCluster::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyCluster) );
      instance.SetNew(&new_TpcPolyCluster);
      instance.SetNewArray(&newArray_TpcPolyCluster);
      instance.SetDelete(&delete_TpcPolyCluster);
      instance.SetDeleteArray(&deleteArray_TpcPolyCluster);
      instance.SetDestructor(&destruct_TpcPolyCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyCluster*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyCluster*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyCluster*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterv1(void *p = nullptr);
   static void *newArray_TpcPolyClusterv1(Long_t size, void *p);
   static void delete_TpcPolyClusterv1(void *p);
   static void deleteArray_TpcPolyClusterv1(void *p);
   static void destruct_TpcPolyClusterv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterv1*)
   {
      ::TpcPolyClusterv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterv1", ::TpcPolyClusterv1::Class_Version(), "TpcPolyClusterv1.h", 8,
                  typeid(::TpcPolyClusterv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterv1) );
      instance.SetNew(&new_TpcPolyClusterv1);
      instance.SetNewArray(&newArray_TpcPolyClusterv1);
      instance.SetDelete(&delete_TpcPolyClusterv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterv1);
      instance.SetDestructor(&destruct_TpcPolyClusterv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterContainer(void *p = nullptr);
   static void *newArray_TpcPolyClusterContainer(Long_t size, void *p);
   static void delete_TpcPolyClusterContainer(void *p);
   static void deleteArray_TpcPolyClusterContainer(void *p);
   static void destruct_TpcPolyClusterContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterContainer*)
   {
      ::TpcPolyClusterContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterContainer", ::TpcPolyClusterContainer::Class_Version(), "TpcPolyClusterContainer.h", 9,
                  typeid(::TpcPolyClusterContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterContainer::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterContainer) );
      instance.SetNew(&new_TpcPolyClusterContainer);
      instance.SetNewArray(&newArray_TpcPolyClusterContainer);
      instance.SetDelete(&delete_TpcPolyClusterContainer);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterContainer);
      instance.SetDestructor(&destruct_TpcPolyClusterContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterContainerv1(void *p = nullptr);
   static void *newArray_TpcPolyClusterContainerv1(Long_t size, void *p);
   static void delete_TpcPolyClusterContainerv1(void *p);
   static void deleteArray_TpcPolyClusterContainerv1(void *p);
   static void destruct_TpcPolyClusterContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterContainerv1*)
   {
      ::TpcPolyClusterContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterContainerv1", ::TpcPolyClusterContainerv1::Class_Version(), "TpcPolyClusterContainerv1.h", 10,
                  typeid(::TpcPolyClusterContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterContainerv1) );
      instance.SetNew(&new_TpcPolyClusterContainerv1);
      instance.SetNewArray(&newArray_TpcPolyClusterContainerv1);
      instance.SetDelete(&delete_TpcPolyClusterContainerv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterContainerv1);
      instance.SetDestructor(&destruct_TpcPolyClusterContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterTrack(void *p = nullptr);
   static void *newArray_TpcPolyClusterTrack(Long_t size, void *p);
   static void delete_TpcPolyClusterTrack(void *p);
   static void deleteArray_TpcPolyClusterTrack(void *p);
   static void destruct_TpcPolyClusterTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterTrack*)
   {
      ::TpcPolyClusterTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterTrack", ::TpcPolyClusterTrack::Class_Version(), "TpcPolyClusterTrack.h", 10,
                  typeid(::TpcPolyClusterTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterTrack::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterTrack) );
      instance.SetNew(&new_TpcPolyClusterTrack);
      instance.SetNewArray(&newArray_TpcPolyClusterTrack);
      instance.SetDelete(&delete_TpcPolyClusterTrack);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterTrack);
      instance.SetDestructor(&destruct_TpcPolyClusterTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterTrackv1(void *p = nullptr);
   static void *newArray_TpcPolyClusterTrackv1(Long_t size, void *p);
   static void delete_TpcPolyClusterTrackv1(void *p);
   static void deleteArray_TpcPolyClusterTrackv1(void *p);
   static void destruct_TpcPolyClusterTrackv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterTrackv1*)
   {
      ::TpcPolyClusterTrackv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterTrackv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterTrackv1", ::TpcPolyClusterTrackv1::Class_Version(), "TpcPolyClusterTrackv1.h", 8,
                  typeid(::TpcPolyClusterTrackv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterTrackv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterTrackv1) );
      instance.SetNew(&new_TpcPolyClusterTrackv1);
      instance.SetNewArray(&newArray_TpcPolyClusterTrackv1);
      instance.SetDelete(&delete_TpcPolyClusterTrackv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterTrackv1);
      instance.SetDestructor(&destruct_TpcPolyClusterTrackv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterTrackv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterTrackv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterTrackv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterTrackContainer(void *p = nullptr);
   static void *newArray_TpcPolyClusterTrackContainer(Long_t size, void *p);
   static void delete_TpcPolyClusterTrackContainer(void *p);
   static void deleteArray_TpcPolyClusterTrackContainer(void *p);
   static void destruct_TpcPolyClusterTrackContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterTrackContainer*)
   {
      ::TpcPolyClusterTrackContainer *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterTrackContainer >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterTrackContainer", ::TpcPolyClusterTrackContainer::Class_Version(), "TpcPolyClusterTrackContainer.h", 9,
                  typeid(::TpcPolyClusterTrackContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterTrackContainer::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterTrackContainer) );
      instance.SetNew(&new_TpcPolyClusterTrackContainer);
      instance.SetNewArray(&newArray_TpcPolyClusterTrackContainer);
      instance.SetDelete(&delete_TpcPolyClusterTrackContainer);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterTrackContainer);
      instance.SetDestructor(&destruct_TpcPolyClusterTrackContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterTrackContainer*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterTrackContainer*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterTrackContainer*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPolyClusterTrackContainerv1(void *p = nullptr);
   static void *newArray_TpcPolyClusterTrackContainerv1(Long_t size, void *p);
   static void delete_TpcPolyClusterTrackContainerv1(void *p);
   static void deleteArray_TpcPolyClusterTrackContainerv1(void *p);
   static void destruct_TpcPolyClusterTrackContainerv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPolyClusterTrackContainerv1*)
   {
      ::TpcPolyClusterTrackContainerv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPolyClusterTrackContainerv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPolyClusterTrackContainerv1", ::TpcPolyClusterTrackContainerv1::Class_Version(), "TpcPolyClusterTrackContainerv1.h", 10,
                  typeid(::TpcPolyClusterTrackContainerv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPolyClusterTrackContainerv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPolyClusterTrackContainerv1) );
      instance.SetNew(&new_TpcPolyClusterTrackContainerv1);
      instance.SetNewArray(&newArray_TpcPolyClusterTrackContainerv1);
      instance.SetDelete(&delete_TpcPolyClusterTrackContainerv1);
      instance.SetDeleteArray(&deleteArray_TpcPolyClusterTrackContainerv1);
      instance.SetDestructor(&destruct_TpcPolyClusterTrackContainerv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPolyClusterTrackContainerv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPolyClusterTrackContainerv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPolyClusterTrackContainerv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *IdealPadMap_Dictionary();
   static void IdealPadMap_TClassManip(TClass*);
   static void *new_IdealPadMap(void *p = nullptr);
   static void *newArray_IdealPadMap(Long_t size, void *p);
   static void delete_IdealPadMap(void *p);
   static void deleteArray_IdealPadMap(void *p);
   static void destruct_IdealPadMap(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::IdealPadMap*)
   {
      ::IdealPadMap *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::IdealPadMap));
      static ::ROOT::TGenericClassInfo 
         instance("IdealPadMap", "IdealPadMap.h", 10,
                  typeid(::IdealPadMap), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &IdealPadMap_Dictionary, isa_proxy, 4,
                  sizeof(::IdealPadMap) );
      instance.SetNew(&new_IdealPadMap);
      instance.SetNewArray(&newArray_IdealPadMap);
      instance.SetDelete(&delete_IdealPadMap);
      instance.SetDeleteArray(&deleteArray_IdealPadMap);
      instance.SetDestructor(&destruct_IdealPadMap);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::IdealPadMap*)
   {
      return GenerateInitInstanceLocal(static_cast<::IdealPadMap*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::IdealPadMap*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *IdealPadMap_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::IdealPadMap*>(nullptr))->GetClass();
      IdealPadMap_TClassManip(theClass);
   return theClass;
   }

   static void IdealPadMap_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPadMap(void *p = nullptr);
   static void *newArray_TpcPadMap(Long_t size, void *p);
   static void delete_TpcPadMap(void *p);
   static void deleteArray_TpcPadMap(void *p);
   static void destruct_TpcPadMap(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPadMap*)
   {
      ::TpcPadMap *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPadMap >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPadMap", ::TpcPadMap::Class_Version(), "TpcPadMap.h", 8,
                  typeid(::TpcPadMap), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPadMap::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPadMap) );
      instance.SetNew(&new_TpcPadMap);
      instance.SetNewArray(&newArray_TpcPadMap);
      instance.SetDelete(&delete_TpcPadMap);
      instance.SetDeleteArray(&deleteArray_TpcPadMap);
      instance.SetDestructor(&destruct_TpcPadMap);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPadMap*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPadMap*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPadMap*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TpcPadMapv1(void *p = nullptr);
   static void *newArray_TpcPadMapv1(Long_t size, void *p);
   static void delete_TpcPadMapv1(void *p);
   static void deleteArray_TpcPadMapv1(void *p);
   static void destruct_TpcPadMapv1(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TpcPadMapv1*)
   {
      ::TpcPadMapv1 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TpcPadMapv1 >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TpcPadMapv1", ::TpcPadMapv1::Class_Version(), "TpcPadMapv1.h", 6,
                  typeid(::TpcPadMapv1), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TpcPadMapv1::Dictionary, isa_proxy, 4,
                  sizeof(::TpcPadMapv1) );
      instance.SetNew(&new_TpcPadMapv1);
      instance.SetNewArray(&newArray_TpcPadMapv1);
      instance.SetDelete(&delete_TpcPadMapv1);
      instance.SetDeleteArray(&deleteArray_TpcPadMapv1);
      instance.SetDestructor(&destruct_TpcPadMapv1);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TpcPadMapv1*)
   {
      return GenerateInitInstanceLocal(static_cast<::TpcPadMapv1*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TpcPadMapv1*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr InModuleTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *InModuleTrack::Class_Name()
{
   return "InModuleTrack";
}

//______________________________________________________________________________
const char *InModuleTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int InModuleTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *InModuleTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *InModuleTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr InModuleTrackv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *InModuleTrackv1::Class_Name()
{
   return "InModuleTrackv1";
}

//______________________________________________________________________________
const char *InModuleTrackv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int InModuleTrackv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *InModuleTrackv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *InModuleTrackv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr InModuleTrackContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *InModuleTrackContainer::Class_Name()
{
   return "InModuleTrackContainer";
}

//______________________________________________________________________________
const char *InModuleTrackContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int InModuleTrackContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *InModuleTrackContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *InModuleTrackContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr InModuleTrackContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *InModuleTrackContainerv1::Class_Name()
{
   return "InModuleTrackContainerv1";
}

//______________________________________________________________________________
const char *InModuleTrackContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int InModuleTrackContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *InModuleTrackContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *InModuleTrackContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InModuleTrackContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrack::Class_Name()
{
   return "FullTrack";
}

//______________________________________________________________________________
const char *FullTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackv1::Class_Name()
{
   return "FullTrackv1";
}

//______________________________________________________________________________
const char *FullTrackv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackContainer::Class_Name()
{
   return "FullTrackContainer";
}

//______________________________________________________________________________
const char *FullTrackContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackContainerv1::Class_Name()
{
   return "FullTrackContainerv1";
}

//______________________________________________________________________________
const char *FullTrackContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackVertex::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackVertex::Class_Name()
{
   return "FullTrackVertex";
}

//______________________________________________________________________________
const char *FullTrackVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertex*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertex*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertex*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackVertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertex*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackVertexv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackVertexv1::Class_Name()
{
   return "FullTrackVertexv1";
}

//______________________________________________________________________________
const char *FullTrackVertexv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackVertexv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackVertexv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackVertexv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackVertexContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackVertexContainer::Class_Name()
{
   return "FullTrackVertexContainer";
}

//______________________________________________________________________________
const char *FullTrackVertexContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackVertexContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackVertexContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackVertexContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FullTrackVertexContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FullTrackVertexContainerv1::Class_Name()
{
   return "FullTrackVertexContainerv1";
}

//______________________________________________________________________________
const char *FullTrackVertexContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FullTrackVertexContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FullTrackVertexContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FullTrackVertexContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FullTrackVertexContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrack::Class_Name()
{
   return "FinalTrack";
}

//______________________________________________________________________________
const char *FinalTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackv1::Class_Name()
{
   return "FinalTrackv1";
}

//______________________________________________________________________________
const char *FinalTrackv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackContainer::Class_Name()
{
   return "FinalTrackContainer";
}

//______________________________________________________________________________
const char *FinalTrackContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackContainerv1::Class_Name()
{
   return "FinalTrackContainerv1";
}

//______________________________________________________________________________
const char *FinalTrackContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackVertex::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackVertex::Class_Name()
{
   return "FinalTrackVertex";
}

//______________________________________________________________________________
const char *FinalTrackVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertex*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertex*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertex*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackVertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertex*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackVertexv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackVertexv1::Class_Name()
{
   return "FinalTrackVertexv1";
}

//______________________________________________________________________________
const char *FinalTrackVertexv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackVertexv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackVertexv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackVertexv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackVertexContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackVertexContainer::Class_Name()
{
   return "FinalTrackVertexContainer";
}

//______________________________________________________________________________
const char *FinalTrackVertexContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackVertexContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackVertexContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackVertexContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FinalTrackVertexContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *FinalTrackVertexContainerv1::Class_Name()
{
   return "FinalTrackVertexContainerv1";
}

//______________________________________________________________________________
const char *FinalTrackVertexContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int FinalTrackVertexContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FinalTrackVertexContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FinalTrackVertexContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FinalTrackVertexContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyTrack::Class_Name()
{
   return "TpcPolyTrack";
}

//______________________________________________________________________________
const char *TpcPolyTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyTrackv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyTrackv1::Class_Name()
{
   return "TpcPolyTrackv1";
}

//______________________________________________________________________________
const char *TpcPolyTrackv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyTrackv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyTrackv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyTrackv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyTrackContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyTrackContainer::Class_Name()
{
   return "TpcPolyTrackContainer";
}

//______________________________________________________________________________
const char *TpcPolyTrackContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyTrackContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyTrackContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyTrackContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyTrackContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyTrackContainerv1::Class_Name()
{
   return "TpcPolyTrackContainerv1";
}

//______________________________________________________________________________
const char *TpcPolyTrackContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyTrackContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyTrackContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyTrackContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyTrackContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyCluster::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyCluster::Class_Name()
{
   return "TpcPolyCluster";
}

//______________________________________________________________________________
const char *TpcPolyCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyCluster*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyCluster*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyCluster*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyCluster::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyCluster*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterv1::Class_Name()
{
   return "TpcPolyClusterv1";
}

//______________________________________________________________________________
const char *TpcPolyClusterv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterContainer::Class_Name()
{
   return "TpcPolyClusterContainer";
}

//______________________________________________________________________________
const char *TpcPolyClusterContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterContainerv1::Class_Name()
{
   return "TpcPolyClusterContainerv1";
}

//______________________________________________________________________________
const char *TpcPolyClusterContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterTrack::Class_Name()
{
   return "TpcPolyClusterTrack";
}

//______________________________________________________________________________
const char *TpcPolyClusterTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterTrackv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterTrackv1::Class_Name()
{
   return "TpcPolyClusterTrackv1";
}

//______________________________________________________________________________
const char *TpcPolyClusterTrackv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterTrackv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterTrackContainer::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterTrackContainer::Class_Name()
{
   return "TpcPolyClusterTrackContainer";
}

//______________________________________________________________________________
const char *TpcPolyClusterTrackContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainer*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterTrackContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainer*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainer*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainer*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPolyClusterTrackContainerv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPolyClusterTrackContainerv1::Class_Name()
{
   return "TpcPolyClusterTrackContainerv1";
}

//______________________________________________________________________________
const char *TpcPolyClusterTrackContainerv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainerv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPolyClusterTrackContainerv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainerv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackContainerv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainerv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPolyClusterTrackContainerv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPolyClusterTrackContainerv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPadMap::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPadMap::Class_Name()
{
   return "TpcPadMap";
}

//______________________________________________________________________________
const char *TpcPadMap::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMap*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPadMap::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMap*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPadMap::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMap*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPadMap::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMap*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TpcPadMapv1::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TpcPadMapv1::Class_Name()
{
   return "TpcPadMapv1";
}

//______________________________________________________________________________
const char *TpcPadMapv1::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMapv1*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TpcPadMapv1::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMapv1*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TpcPadMapv1::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMapv1*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TpcPadMapv1::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TpcPadMapv1*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void InModuleTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class InModuleTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(InModuleTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(InModuleTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_InModuleTrack(void *p) {
      return  p ? new(p) ::InModuleTrack : new ::InModuleTrack;
   }
   static void *newArray_InModuleTrack(Long_t nElements, void *p) {
      return p ? new(p) ::InModuleTrack[nElements] : new ::InModuleTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_InModuleTrack(void *p) {
      delete (static_cast<::InModuleTrack*>(p));
   }
   static void deleteArray_InModuleTrack(void *p) {
      delete [] (static_cast<::InModuleTrack*>(p));
   }
   static void destruct_InModuleTrack(void *p) {
      typedef ::InModuleTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::InModuleTrack

//______________________________________________________________________________
void InModuleTrackv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class InModuleTrackv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(InModuleTrackv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(InModuleTrackv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_InModuleTrackv1(void *p) {
      return  p ? new(p) ::InModuleTrackv1 : new ::InModuleTrackv1;
   }
   static void *newArray_InModuleTrackv1(Long_t nElements, void *p) {
      return p ? new(p) ::InModuleTrackv1[nElements] : new ::InModuleTrackv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_InModuleTrackv1(void *p) {
      delete (static_cast<::InModuleTrackv1*>(p));
   }
   static void deleteArray_InModuleTrackv1(void *p) {
      delete [] (static_cast<::InModuleTrackv1*>(p));
   }
   static void destruct_InModuleTrackv1(void *p) {
      typedef ::InModuleTrackv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::InModuleTrackv1

//______________________________________________________________________________
void InModuleTrackContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class InModuleTrackContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(InModuleTrackContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(InModuleTrackContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_InModuleTrackContainer(void *p) {
      return  p ? new(p) ::InModuleTrackContainer : new ::InModuleTrackContainer;
   }
   static void *newArray_InModuleTrackContainer(Long_t nElements, void *p) {
      return p ? new(p) ::InModuleTrackContainer[nElements] : new ::InModuleTrackContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_InModuleTrackContainer(void *p) {
      delete (static_cast<::InModuleTrackContainer*>(p));
   }
   static void deleteArray_InModuleTrackContainer(void *p) {
      delete [] (static_cast<::InModuleTrackContainer*>(p));
   }
   static void destruct_InModuleTrackContainer(void *p) {
      typedef ::InModuleTrackContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::InModuleTrackContainer

//______________________________________________________________________________
void InModuleTrackContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class InModuleTrackContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(InModuleTrackContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(InModuleTrackContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_InModuleTrackContainerv1(void *p) {
      return  p ? new(p) ::InModuleTrackContainerv1 : new ::InModuleTrackContainerv1;
   }
   static void *newArray_InModuleTrackContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::InModuleTrackContainerv1[nElements] : new ::InModuleTrackContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_InModuleTrackContainerv1(void *p) {
      delete (static_cast<::InModuleTrackContainerv1*>(p));
   }
   static void deleteArray_InModuleTrackContainerv1(void *p) {
      delete [] (static_cast<::InModuleTrackContainerv1*>(p));
   }
   static void destruct_InModuleTrackContainerv1(void *p) {
      typedef ::InModuleTrackContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::InModuleTrackContainerv1

//______________________________________________________________________________
void FullTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrack(void *p) {
      return  p ? new(p) ::FullTrack : new ::FullTrack;
   }
   static void *newArray_FullTrack(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrack[nElements] : new ::FullTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrack(void *p) {
      delete (static_cast<::FullTrack*>(p));
   }
   static void deleteArray_FullTrack(void *p) {
      delete [] (static_cast<::FullTrack*>(p));
   }
   static void destruct_FullTrack(void *p) {
      typedef ::FullTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrack

//______________________________________________________________________________
void FullTrackv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackv1(void *p) {
      return  p ? new(p) ::FullTrackv1 : new ::FullTrackv1;
   }
   static void *newArray_FullTrackv1(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackv1[nElements] : new ::FullTrackv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackv1(void *p) {
      delete (static_cast<::FullTrackv1*>(p));
   }
   static void deleteArray_FullTrackv1(void *p) {
      delete [] (static_cast<::FullTrackv1*>(p));
   }
   static void destruct_FullTrackv1(void *p) {
      typedef ::FullTrackv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackv1

//______________________________________________________________________________
void FullTrackContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackContainer(void *p) {
      return  p ? new(p) ::FullTrackContainer : new ::FullTrackContainer;
   }
   static void *newArray_FullTrackContainer(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackContainer[nElements] : new ::FullTrackContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackContainer(void *p) {
      delete (static_cast<::FullTrackContainer*>(p));
   }
   static void deleteArray_FullTrackContainer(void *p) {
      delete [] (static_cast<::FullTrackContainer*>(p));
   }
   static void destruct_FullTrackContainer(void *p) {
      typedef ::FullTrackContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackContainer

//______________________________________________________________________________
void FullTrackContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackContainerv1(void *p) {
      return  p ? new(p) ::FullTrackContainerv1 : new ::FullTrackContainerv1;
   }
   static void *newArray_FullTrackContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackContainerv1[nElements] : new ::FullTrackContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackContainerv1(void *p) {
      delete (static_cast<::FullTrackContainerv1*>(p));
   }
   static void deleteArray_FullTrackContainerv1(void *p) {
      delete [] (static_cast<::FullTrackContainerv1*>(p));
   }
   static void destruct_FullTrackContainerv1(void *p) {
      typedef ::FullTrackContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackContainerv1

//______________________________________________________________________________
void FullTrackVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackVertex::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackVertex(void *p) {
      return  p ? new(p) ::FullTrackVertex : new ::FullTrackVertex;
   }
   static void *newArray_FullTrackVertex(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackVertex[nElements] : new ::FullTrackVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackVertex(void *p) {
      delete (static_cast<::FullTrackVertex*>(p));
   }
   static void deleteArray_FullTrackVertex(void *p) {
      delete [] (static_cast<::FullTrackVertex*>(p));
   }
   static void destruct_FullTrackVertex(void *p) {
      typedef ::FullTrackVertex current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackVertex

//______________________________________________________________________________
void FullTrackVertexv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackVertexv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackVertexv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackVertexv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackVertexv1(void *p) {
      return  p ? new(p) ::FullTrackVertexv1 : new ::FullTrackVertexv1;
   }
   static void *newArray_FullTrackVertexv1(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackVertexv1[nElements] : new ::FullTrackVertexv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackVertexv1(void *p) {
      delete (static_cast<::FullTrackVertexv1*>(p));
   }
   static void deleteArray_FullTrackVertexv1(void *p) {
      delete [] (static_cast<::FullTrackVertexv1*>(p));
   }
   static void destruct_FullTrackVertexv1(void *p) {
      typedef ::FullTrackVertexv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackVertexv1

//______________________________________________________________________________
void FullTrackVertexContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackVertexContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackVertexContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackVertexContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackVertexContainer(void *p) {
      return  p ? new(p) ::FullTrackVertexContainer : new ::FullTrackVertexContainer;
   }
   static void *newArray_FullTrackVertexContainer(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackVertexContainer[nElements] : new ::FullTrackVertexContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackVertexContainer(void *p) {
      delete (static_cast<::FullTrackVertexContainer*>(p));
   }
   static void deleteArray_FullTrackVertexContainer(void *p) {
      delete [] (static_cast<::FullTrackVertexContainer*>(p));
   }
   static void destruct_FullTrackVertexContainer(void *p) {
      typedef ::FullTrackVertexContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackVertexContainer

//______________________________________________________________________________
void FullTrackVertexContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FullTrackVertexContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FullTrackVertexContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FullTrackVertexContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FullTrackVertexContainerv1(void *p) {
      return  p ? new(p) ::FullTrackVertexContainerv1 : new ::FullTrackVertexContainerv1;
   }
   static void *newArray_FullTrackVertexContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::FullTrackVertexContainerv1[nElements] : new ::FullTrackVertexContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FullTrackVertexContainerv1(void *p) {
      delete (static_cast<::FullTrackVertexContainerv1*>(p));
   }
   static void deleteArray_FullTrackVertexContainerv1(void *p) {
      delete [] (static_cast<::FullTrackVertexContainerv1*>(p));
   }
   static void destruct_FullTrackVertexContainerv1(void *p) {
      typedef ::FullTrackVertexContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FullTrackVertexContainerv1

//______________________________________________________________________________
void FinalTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrack(void *p) {
      return  p ? new(p) ::FinalTrack : new ::FinalTrack;
   }
   static void *newArray_FinalTrack(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrack[nElements] : new ::FinalTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrack(void *p) {
      delete (static_cast<::FinalTrack*>(p));
   }
   static void deleteArray_FinalTrack(void *p) {
      delete [] (static_cast<::FinalTrack*>(p));
   }
   static void destruct_FinalTrack(void *p) {
      typedef ::FinalTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrack

//______________________________________________________________________________
void FinalTrackv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackv1(void *p) {
      return  p ? new(p) ::FinalTrackv1 : new ::FinalTrackv1;
   }
   static void *newArray_FinalTrackv1(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackv1[nElements] : new ::FinalTrackv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackv1(void *p) {
      delete (static_cast<::FinalTrackv1*>(p));
   }
   static void deleteArray_FinalTrackv1(void *p) {
      delete [] (static_cast<::FinalTrackv1*>(p));
   }
   static void destruct_FinalTrackv1(void *p) {
      typedef ::FinalTrackv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackv1

//______________________________________________________________________________
void FinalTrackContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackContainer(void *p) {
      return  p ? new(p) ::FinalTrackContainer : new ::FinalTrackContainer;
   }
   static void *newArray_FinalTrackContainer(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackContainer[nElements] : new ::FinalTrackContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackContainer(void *p) {
      delete (static_cast<::FinalTrackContainer*>(p));
   }
   static void deleteArray_FinalTrackContainer(void *p) {
      delete [] (static_cast<::FinalTrackContainer*>(p));
   }
   static void destruct_FinalTrackContainer(void *p) {
      typedef ::FinalTrackContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackContainer

//______________________________________________________________________________
void FinalTrackContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackContainerv1(void *p) {
      return  p ? new(p) ::FinalTrackContainerv1 : new ::FinalTrackContainerv1;
   }
   static void *newArray_FinalTrackContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackContainerv1[nElements] : new ::FinalTrackContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackContainerv1(void *p) {
      delete (static_cast<::FinalTrackContainerv1*>(p));
   }
   static void deleteArray_FinalTrackContainerv1(void *p) {
      delete [] (static_cast<::FinalTrackContainerv1*>(p));
   }
   static void destruct_FinalTrackContainerv1(void *p) {
      typedef ::FinalTrackContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackContainerv1

//______________________________________________________________________________
void FinalTrackVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackVertex::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackVertex(void *p) {
      return  p ? new(p) ::FinalTrackVertex : new ::FinalTrackVertex;
   }
   static void *newArray_FinalTrackVertex(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackVertex[nElements] : new ::FinalTrackVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackVertex(void *p) {
      delete (static_cast<::FinalTrackVertex*>(p));
   }
   static void deleteArray_FinalTrackVertex(void *p) {
      delete [] (static_cast<::FinalTrackVertex*>(p));
   }
   static void destruct_FinalTrackVertex(void *p) {
      typedef ::FinalTrackVertex current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackVertex

//______________________________________________________________________________
void FinalTrackVertexv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackVertexv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackVertexv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackVertexv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackVertexv1(void *p) {
      return  p ? new(p) ::FinalTrackVertexv1 : new ::FinalTrackVertexv1;
   }
   static void *newArray_FinalTrackVertexv1(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackVertexv1[nElements] : new ::FinalTrackVertexv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackVertexv1(void *p) {
      delete (static_cast<::FinalTrackVertexv1*>(p));
   }
   static void deleteArray_FinalTrackVertexv1(void *p) {
      delete [] (static_cast<::FinalTrackVertexv1*>(p));
   }
   static void destruct_FinalTrackVertexv1(void *p) {
      typedef ::FinalTrackVertexv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackVertexv1

//______________________________________________________________________________
void FinalTrackVertexContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackVertexContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackVertexContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackVertexContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackVertexContainer(void *p) {
      return  p ? new(p) ::FinalTrackVertexContainer : new ::FinalTrackVertexContainer;
   }
   static void *newArray_FinalTrackVertexContainer(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackVertexContainer[nElements] : new ::FinalTrackVertexContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackVertexContainer(void *p) {
      delete (static_cast<::FinalTrackVertexContainer*>(p));
   }
   static void deleteArray_FinalTrackVertexContainer(void *p) {
      delete [] (static_cast<::FinalTrackVertexContainer*>(p));
   }
   static void destruct_FinalTrackVertexContainer(void *p) {
      typedef ::FinalTrackVertexContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackVertexContainer

//______________________________________________________________________________
void FinalTrackVertexContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class FinalTrackVertexContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FinalTrackVertexContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(FinalTrackVertexContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FinalTrackVertexContainerv1(void *p) {
      return  p ? new(p) ::FinalTrackVertexContainerv1 : new ::FinalTrackVertexContainerv1;
   }
   static void *newArray_FinalTrackVertexContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::FinalTrackVertexContainerv1[nElements] : new ::FinalTrackVertexContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_FinalTrackVertexContainerv1(void *p) {
      delete (static_cast<::FinalTrackVertexContainerv1*>(p));
   }
   static void deleteArray_FinalTrackVertexContainerv1(void *p) {
      delete [] (static_cast<::FinalTrackVertexContainerv1*>(p));
   }
   static void destruct_FinalTrackVertexContainerv1(void *p) {
      typedef ::FinalTrackVertexContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FinalTrackVertexContainerv1

//______________________________________________________________________________
void TpcPolyTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyTrack(void *p) {
      return  p ? new(p) ::TpcPolyTrack : new ::TpcPolyTrack;
   }
   static void *newArray_TpcPolyTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyTrack[nElements] : new ::TpcPolyTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyTrack(void *p) {
      delete (static_cast<::TpcPolyTrack*>(p));
   }
   static void deleteArray_TpcPolyTrack(void *p) {
      delete [] (static_cast<::TpcPolyTrack*>(p));
   }
   static void destruct_TpcPolyTrack(void *p) {
      typedef ::TpcPolyTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyTrack

//______________________________________________________________________________
void TpcPolyTrackv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyTrackv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyTrackv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyTrackv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyTrackv1(void *p) {
      return  p ? new(p) ::TpcPolyTrackv1 : new ::TpcPolyTrackv1;
   }
   static void *newArray_TpcPolyTrackv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyTrackv1[nElements] : new ::TpcPolyTrackv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyTrackv1(void *p) {
      delete (static_cast<::TpcPolyTrackv1*>(p));
   }
   static void deleteArray_TpcPolyTrackv1(void *p) {
      delete [] (static_cast<::TpcPolyTrackv1*>(p));
   }
   static void destruct_TpcPolyTrackv1(void *p) {
      typedef ::TpcPolyTrackv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyTrackv1

//______________________________________________________________________________
void TpcPolyTrackContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyTrackContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyTrackContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyTrackContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyTrackContainer(void *p) {
      return  p ? new(p) ::TpcPolyTrackContainer : new ::TpcPolyTrackContainer;
   }
   static void *newArray_TpcPolyTrackContainer(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyTrackContainer[nElements] : new ::TpcPolyTrackContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyTrackContainer(void *p) {
      delete (static_cast<::TpcPolyTrackContainer*>(p));
   }
   static void deleteArray_TpcPolyTrackContainer(void *p) {
      delete [] (static_cast<::TpcPolyTrackContainer*>(p));
   }
   static void destruct_TpcPolyTrackContainer(void *p) {
      typedef ::TpcPolyTrackContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyTrackContainer

//______________________________________________________________________________
void TpcPolyTrackContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyTrackContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyTrackContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyTrackContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyTrackContainerv1(void *p) {
      return  p ? new(p) ::TpcPolyTrackContainerv1 : new ::TpcPolyTrackContainerv1;
   }
   static void *newArray_TpcPolyTrackContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyTrackContainerv1[nElements] : new ::TpcPolyTrackContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyTrackContainerv1(void *p) {
      delete (static_cast<::TpcPolyTrackContainerv1*>(p));
   }
   static void deleteArray_TpcPolyTrackContainerv1(void *p) {
      delete [] (static_cast<::TpcPolyTrackContainerv1*>(p));
   }
   static void destruct_TpcPolyTrackContainerv1(void *p) {
      typedef ::TpcPolyTrackContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyTrackContainerv1

//______________________________________________________________________________
void TpcPolyCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyCluster::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyCluster(void *p) {
      return  p ? new(p) ::TpcPolyCluster : new ::TpcPolyCluster;
   }
   static void *newArray_TpcPolyCluster(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyCluster[nElements] : new ::TpcPolyCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyCluster(void *p) {
      delete (static_cast<::TpcPolyCluster*>(p));
   }
   static void deleteArray_TpcPolyCluster(void *p) {
      delete [] (static_cast<::TpcPolyCluster*>(p));
   }
   static void destruct_TpcPolyCluster(void *p) {
      typedef ::TpcPolyCluster current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyCluster

//______________________________________________________________________________
void TpcPolyClusterv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterv1(void *p) {
      return  p ? new(p) ::TpcPolyClusterv1 : new ::TpcPolyClusterv1;
   }
   static void *newArray_TpcPolyClusterv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterv1[nElements] : new ::TpcPolyClusterv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterv1(void *p) {
      delete (static_cast<::TpcPolyClusterv1*>(p));
   }
   static void deleteArray_TpcPolyClusterv1(void *p) {
      delete [] (static_cast<::TpcPolyClusterv1*>(p));
   }
   static void destruct_TpcPolyClusterv1(void *p) {
      typedef ::TpcPolyClusterv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterv1

//______________________________________________________________________________
void TpcPolyClusterContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterContainer(void *p) {
      return  p ? new(p) ::TpcPolyClusterContainer : new ::TpcPolyClusterContainer;
   }
   static void *newArray_TpcPolyClusterContainer(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterContainer[nElements] : new ::TpcPolyClusterContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterContainer(void *p) {
      delete (static_cast<::TpcPolyClusterContainer*>(p));
   }
   static void deleteArray_TpcPolyClusterContainer(void *p) {
      delete [] (static_cast<::TpcPolyClusterContainer*>(p));
   }
   static void destruct_TpcPolyClusterContainer(void *p) {
      typedef ::TpcPolyClusterContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterContainer

//______________________________________________________________________________
void TpcPolyClusterContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterContainerv1(void *p) {
      return  p ? new(p) ::TpcPolyClusterContainerv1 : new ::TpcPolyClusterContainerv1;
   }
   static void *newArray_TpcPolyClusterContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterContainerv1[nElements] : new ::TpcPolyClusterContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterContainerv1(void *p) {
      delete (static_cast<::TpcPolyClusterContainerv1*>(p));
   }
   static void deleteArray_TpcPolyClusterContainerv1(void *p) {
      delete [] (static_cast<::TpcPolyClusterContainerv1*>(p));
   }
   static void destruct_TpcPolyClusterContainerv1(void *p) {
      typedef ::TpcPolyClusterContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterContainerv1

//______________________________________________________________________________
void TpcPolyClusterTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterTrack(void *p) {
      return  p ? new(p) ::TpcPolyClusterTrack : new ::TpcPolyClusterTrack;
   }
   static void *newArray_TpcPolyClusterTrack(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterTrack[nElements] : new ::TpcPolyClusterTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterTrack(void *p) {
      delete (static_cast<::TpcPolyClusterTrack*>(p));
   }
   static void deleteArray_TpcPolyClusterTrack(void *p) {
      delete [] (static_cast<::TpcPolyClusterTrack*>(p));
   }
   static void destruct_TpcPolyClusterTrack(void *p) {
      typedef ::TpcPolyClusterTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterTrack

//______________________________________________________________________________
void TpcPolyClusterTrackv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterTrackv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterTrackv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterTrackv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterTrackv1(void *p) {
      return  p ? new(p) ::TpcPolyClusterTrackv1 : new ::TpcPolyClusterTrackv1;
   }
   static void *newArray_TpcPolyClusterTrackv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterTrackv1[nElements] : new ::TpcPolyClusterTrackv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterTrackv1(void *p) {
      delete (static_cast<::TpcPolyClusterTrackv1*>(p));
   }
   static void deleteArray_TpcPolyClusterTrackv1(void *p) {
      delete [] (static_cast<::TpcPolyClusterTrackv1*>(p));
   }
   static void destruct_TpcPolyClusterTrackv1(void *p) {
      typedef ::TpcPolyClusterTrackv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterTrackv1

//______________________________________________________________________________
void TpcPolyClusterTrackContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterTrackContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterTrackContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterTrackContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterTrackContainer(void *p) {
      return  p ? new(p) ::TpcPolyClusterTrackContainer : new ::TpcPolyClusterTrackContainer;
   }
   static void *newArray_TpcPolyClusterTrackContainer(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterTrackContainer[nElements] : new ::TpcPolyClusterTrackContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterTrackContainer(void *p) {
      delete (static_cast<::TpcPolyClusterTrackContainer*>(p));
   }
   static void deleteArray_TpcPolyClusterTrackContainer(void *p) {
      delete [] (static_cast<::TpcPolyClusterTrackContainer*>(p));
   }
   static void destruct_TpcPolyClusterTrackContainer(void *p) {
      typedef ::TpcPolyClusterTrackContainer current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterTrackContainer

//______________________________________________________________________________
void TpcPolyClusterTrackContainerv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPolyClusterTrackContainerv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPolyClusterTrackContainerv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPolyClusterTrackContainerv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPolyClusterTrackContainerv1(void *p) {
      return  p ? new(p) ::TpcPolyClusterTrackContainerv1 : new ::TpcPolyClusterTrackContainerv1;
   }
   static void *newArray_TpcPolyClusterTrackContainerv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPolyClusterTrackContainerv1[nElements] : new ::TpcPolyClusterTrackContainerv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPolyClusterTrackContainerv1(void *p) {
      delete (static_cast<::TpcPolyClusterTrackContainerv1*>(p));
   }
   static void deleteArray_TpcPolyClusterTrackContainerv1(void *p) {
      delete [] (static_cast<::TpcPolyClusterTrackContainerv1*>(p));
   }
   static void destruct_TpcPolyClusterTrackContainerv1(void *p) {
      typedef ::TpcPolyClusterTrackContainerv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPolyClusterTrackContainerv1

namespace ROOT {
   // Wrappers around operator new
   static void *new_IdealPadMap(void *p) {
      return  p ? new(p) ::IdealPadMap : new ::IdealPadMap;
   }
   static void *newArray_IdealPadMap(Long_t nElements, void *p) {
      return p ? new(p) ::IdealPadMap[nElements] : new ::IdealPadMap[nElements];
   }
   // Wrapper around operator delete
   static void delete_IdealPadMap(void *p) {
      delete (static_cast<::IdealPadMap*>(p));
   }
   static void deleteArray_IdealPadMap(void *p) {
      delete [] (static_cast<::IdealPadMap*>(p));
   }
   static void destruct_IdealPadMap(void *p) {
      typedef ::IdealPadMap current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::IdealPadMap

//______________________________________________________________________________
void TpcPadMap::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPadMap.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPadMap::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPadMap::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPadMap(void *p) {
      return  p ? new(p) ::TpcPadMap : new ::TpcPadMap;
   }
   static void *newArray_TpcPadMap(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPadMap[nElements] : new ::TpcPadMap[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPadMap(void *p) {
      delete (static_cast<::TpcPadMap*>(p));
   }
   static void deleteArray_TpcPadMap(void *p) {
      delete [] (static_cast<::TpcPadMap*>(p));
   }
   static void destruct_TpcPadMap(void *p) {
      typedef ::TpcPadMap current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPadMap

//______________________________________________________________________________
void TpcPadMapv1::Streamer(TBuffer &R__b)
{
   // Stream an object of class TpcPadMapv1.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TpcPadMapv1::Class(),this);
   } else {
      R__b.WriteClassBuffer(TpcPadMapv1::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TpcPadMapv1(void *p) {
      return  p ? new(p) ::TpcPadMapv1 : new ::TpcPadMapv1;
   }
   static void *newArray_TpcPadMapv1(Long_t nElements, void *p) {
      return p ? new(p) ::TpcPadMapv1[nElements] : new ::TpcPadMapv1[nElements];
   }
   // Wrapper around operator delete
   static void delete_TpcPadMapv1(void *p) {
      delete (static_cast<::TpcPadMapv1*>(p));
   }
   static void deleteArray_TpcPadMapv1(void *p) {
      delete [] (static_cast<::TpcPadMapv1*>(p));
   }
   static void destruct_TpcPadMapv1(void *p) {
      typedef ::TpcPadMapv1 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TpcPadMapv1

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "vector", 428,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned int>","std::vector<unsigned int, std::allocator<unsigned int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned int>*>(nullptr))->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete (static_cast<vector<unsigned int>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] (static_cast<vector<unsigned int>*>(p));
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary();
   static void vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);
   static void destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<unsigned int,unsigned int> >*)
   {
      vector<pair<unsigned int,unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<unsigned int,unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<unsigned int,unsigned int> >", -2, "vector", 428,
                  typeid(vector<pair<unsigned int,unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<pair<unsigned int,unsigned int> >) );
      instance.SetNew(&new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<unsigned int,unsigned int> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<pair<unsigned int,unsigned int> >","std::vector<std::pair<unsigned int, unsigned int>, std::allocator<std::pair<unsigned int, unsigned int> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<pair<unsigned int,unsigned int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<pair<unsigned int,unsigned int> >*>(nullptr))->GetClass();
      vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<unsigned int,unsigned int> > : new vector<pair<unsigned int,unsigned int> >;
   }
   static void *newArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<unsigned int,unsigned int> >[nElements] : new vector<pair<unsigned int,unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      delete (static_cast<vector<pair<unsigned int,unsigned int> >*>(p));
   }
   static void deleteArray_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      delete [] (static_cast<vector<pair<unsigned int,unsigned int> >*>(p));
   }
   static void destruct_vectorlEpairlEunsignedsPintcOunsignedsPintgRsPgR(void *p) {
      typedef vector<pair<unsigned int,unsigned int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<pair<unsigned int,unsigned int> >

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 428,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETpcPolyTrackmUgR_Dictionary();
   static void vectorlETpcPolyTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlETpcPolyTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlETpcPolyTrackmUgR(Long_t size, void *p);
   static void delete_vectorlETpcPolyTrackmUgR(void *p);
   static void deleteArray_vectorlETpcPolyTrackmUgR(void *p);
   static void destruct_vectorlETpcPolyTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TpcPolyTrack*>*)
   {
      vector<TpcPolyTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TpcPolyTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TpcPolyTrack*>", -2, "vector", 428,
                  typeid(vector<TpcPolyTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETpcPolyTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TpcPolyTrack*>) );
      instance.SetNew(&new_vectorlETpcPolyTrackmUgR);
      instance.SetNewArray(&newArray_vectorlETpcPolyTrackmUgR);
      instance.SetDelete(&delete_vectorlETpcPolyTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETpcPolyTrackmUgR);
      instance.SetDestructor(&destruct_vectorlETpcPolyTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TpcPolyTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TpcPolyTrack*>","std::vector<TpcPolyTrack*, std::allocator<TpcPolyTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TpcPolyTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETpcPolyTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TpcPolyTrack*>*>(nullptr))->GetClass();
      vectorlETpcPolyTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETpcPolyTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETpcPolyTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyTrack*> : new vector<TpcPolyTrack*>;
   }
   static void *newArray_vectorlETpcPolyTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyTrack*>[nElements] : new vector<TpcPolyTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETpcPolyTrackmUgR(void *p) {
      delete (static_cast<vector<TpcPolyTrack*>*>(p));
   }
   static void deleteArray_vectorlETpcPolyTrackmUgR(void *p) {
      delete [] (static_cast<vector<TpcPolyTrack*>*>(p));
   }
   static void destruct_vectorlETpcPolyTrackmUgR(void *p) {
      typedef vector<TpcPolyTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TpcPolyTrack*>

namespace ROOT {
   static TClass *vectorlETpcPolyClusterTrackmUgR_Dictionary();
   static void vectorlETpcPolyClusterTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlETpcPolyClusterTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlETpcPolyClusterTrackmUgR(Long_t size, void *p);
   static void delete_vectorlETpcPolyClusterTrackmUgR(void *p);
   static void deleteArray_vectorlETpcPolyClusterTrackmUgR(void *p);
   static void destruct_vectorlETpcPolyClusterTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TpcPolyClusterTrack*>*)
   {
      vector<TpcPolyClusterTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TpcPolyClusterTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TpcPolyClusterTrack*>", -2, "vector", 428,
                  typeid(vector<TpcPolyClusterTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETpcPolyClusterTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TpcPolyClusterTrack*>) );
      instance.SetNew(&new_vectorlETpcPolyClusterTrackmUgR);
      instance.SetNewArray(&newArray_vectorlETpcPolyClusterTrackmUgR);
      instance.SetDelete(&delete_vectorlETpcPolyClusterTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETpcPolyClusterTrackmUgR);
      instance.SetDestructor(&destruct_vectorlETpcPolyClusterTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TpcPolyClusterTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TpcPolyClusterTrack*>","std::vector<TpcPolyClusterTrack*, std::allocator<TpcPolyClusterTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TpcPolyClusterTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETpcPolyClusterTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TpcPolyClusterTrack*>*>(nullptr))->GetClass();
      vectorlETpcPolyClusterTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETpcPolyClusterTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETpcPolyClusterTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyClusterTrack*> : new vector<TpcPolyClusterTrack*>;
   }
   static void *newArray_vectorlETpcPolyClusterTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyClusterTrack*>[nElements] : new vector<TpcPolyClusterTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETpcPolyClusterTrackmUgR(void *p) {
      delete (static_cast<vector<TpcPolyClusterTrack*>*>(p));
   }
   static void deleteArray_vectorlETpcPolyClusterTrackmUgR(void *p) {
      delete [] (static_cast<vector<TpcPolyClusterTrack*>*>(p));
   }
   static void destruct_vectorlETpcPolyClusterTrackmUgR(void *p) {
      typedef vector<TpcPolyClusterTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TpcPolyClusterTrack*>

namespace ROOT {
   static TClass *vectorlETpcPolyClustermUgR_Dictionary();
   static void vectorlETpcPolyClustermUgR_TClassManip(TClass*);
   static void *new_vectorlETpcPolyClustermUgR(void *p = nullptr);
   static void *newArray_vectorlETpcPolyClustermUgR(Long_t size, void *p);
   static void delete_vectorlETpcPolyClustermUgR(void *p);
   static void deleteArray_vectorlETpcPolyClustermUgR(void *p);
   static void destruct_vectorlETpcPolyClustermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TpcPolyCluster*>*)
   {
      vector<TpcPolyCluster*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TpcPolyCluster*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TpcPolyCluster*>", -2, "vector", 428,
                  typeid(vector<TpcPolyCluster*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETpcPolyClustermUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TpcPolyCluster*>) );
      instance.SetNew(&new_vectorlETpcPolyClustermUgR);
      instance.SetNewArray(&newArray_vectorlETpcPolyClustermUgR);
      instance.SetDelete(&delete_vectorlETpcPolyClustermUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETpcPolyClustermUgR);
      instance.SetDestructor(&destruct_vectorlETpcPolyClustermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TpcPolyCluster*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TpcPolyCluster*>","std::vector<TpcPolyCluster*, std::allocator<TpcPolyCluster*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TpcPolyCluster*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETpcPolyClustermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TpcPolyCluster*>*>(nullptr))->GetClass();
      vectorlETpcPolyClustermUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETpcPolyClustermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETpcPolyClustermUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyCluster*> : new vector<TpcPolyCluster*>;
   }
   static void *newArray_vectorlETpcPolyClustermUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TpcPolyCluster*>[nElements] : new vector<TpcPolyCluster*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETpcPolyClustermUgR(void *p) {
      delete (static_cast<vector<TpcPolyCluster*>*>(p));
   }
   static void deleteArray_vectorlETpcPolyClustermUgR(void *p) {
      delete [] (static_cast<vector<TpcPolyCluster*>*>(p));
   }
   static void destruct_vectorlETpcPolyClustermUgR(void *p) {
      typedef vector<TpcPolyCluster*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TpcPolyCluster*>

namespace ROOT {
   static TClass *vectorlEInModuleTrackmUgR_Dictionary();
   static void vectorlEInModuleTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlEInModuleTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlEInModuleTrackmUgR(Long_t size, void *p);
   static void delete_vectorlEInModuleTrackmUgR(void *p);
   static void deleteArray_vectorlEInModuleTrackmUgR(void *p);
   static void destruct_vectorlEInModuleTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<InModuleTrack*>*)
   {
      vector<InModuleTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<InModuleTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<InModuleTrack*>", -2, "vector", 428,
                  typeid(vector<InModuleTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEInModuleTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<InModuleTrack*>) );
      instance.SetNew(&new_vectorlEInModuleTrackmUgR);
      instance.SetNewArray(&newArray_vectorlEInModuleTrackmUgR);
      instance.SetDelete(&delete_vectorlEInModuleTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEInModuleTrackmUgR);
      instance.SetDestructor(&destruct_vectorlEInModuleTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<InModuleTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<InModuleTrack*>","std::vector<InModuleTrack*, std::allocator<InModuleTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<InModuleTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEInModuleTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<InModuleTrack*>*>(nullptr))->GetClass();
      vectorlEInModuleTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEInModuleTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEInModuleTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<InModuleTrack*> : new vector<InModuleTrack*>;
   }
   static void *newArray_vectorlEInModuleTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<InModuleTrack*>[nElements] : new vector<InModuleTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEInModuleTrackmUgR(void *p) {
      delete (static_cast<vector<InModuleTrack*>*>(p));
   }
   static void deleteArray_vectorlEInModuleTrackmUgR(void *p) {
      delete [] (static_cast<vector<InModuleTrack*>*>(p));
   }
   static void destruct_vectorlEInModuleTrackmUgR(void *p) {
      typedef vector<InModuleTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<InModuleTrack*>

namespace ROOT {
   static TClass *vectorlEFullTrackVertexmUgR_Dictionary();
   static void vectorlEFullTrackVertexmUgR_TClassManip(TClass*);
   static void *new_vectorlEFullTrackVertexmUgR(void *p = nullptr);
   static void *newArray_vectorlEFullTrackVertexmUgR(Long_t size, void *p);
   static void delete_vectorlEFullTrackVertexmUgR(void *p);
   static void deleteArray_vectorlEFullTrackVertexmUgR(void *p);
   static void destruct_vectorlEFullTrackVertexmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<FullTrackVertex*>*)
   {
      vector<FullTrackVertex*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<FullTrackVertex*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<FullTrackVertex*>", -2, "vector", 428,
                  typeid(vector<FullTrackVertex*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEFullTrackVertexmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<FullTrackVertex*>) );
      instance.SetNew(&new_vectorlEFullTrackVertexmUgR);
      instance.SetNewArray(&newArray_vectorlEFullTrackVertexmUgR);
      instance.SetDelete(&delete_vectorlEFullTrackVertexmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEFullTrackVertexmUgR);
      instance.SetDestructor(&destruct_vectorlEFullTrackVertexmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<FullTrackVertex*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<FullTrackVertex*>","std::vector<FullTrackVertex*, std::allocator<FullTrackVertex*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<FullTrackVertex*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEFullTrackVertexmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<FullTrackVertex*>*>(nullptr))->GetClass();
      vectorlEFullTrackVertexmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEFullTrackVertexmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEFullTrackVertexmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FullTrackVertex*> : new vector<FullTrackVertex*>;
   }
   static void *newArray_vectorlEFullTrackVertexmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FullTrackVertex*>[nElements] : new vector<FullTrackVertex*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEFullTrackVertexmUgR(void *p) {
      delete (static_cast<vector<FullTrackVertex*>*>(p));
   }
   static void deleteArray_vectorlEFullTrackVertexmUgR(void *p) {
      delete [] (static_cast<vector<FullTrackVertex*>*>(p));
   }
   static void destruct_vectorlEFullTrackVertexmUgR(void *p) {
      typedef vector<FullTrackVertex*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<FullTrackVertex*>

namespace ROOT {
   static TClass *vectorlEFullTrackmUgR_Dictionary();
   static void vectorlEFullTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlEFullTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlEFullTrackmUgR(Long_t size, void *p);
   static void delete_vectorlEFullTrackmUgR(void *p);
   static void deleteArray_vectorlEFullTrackmUgR(void *p);
   static void destruct_vectorlEFullTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<FullTrack*>*)
   {
      vector<FullTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<FullTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<FullTrack*>", -2, "vector", 428,
                  typeid(vector<FullTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEFullTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<FullTrack*>) );
      instance.SetNew(&new_vectorlEFullTrackmUgR);
      instance.SetNewArray(&newArray_vectorlEFullTrackmUgR);
      instance.SetDelete(&delete_vectorlEFullTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEFullTrackmUgR);
      instance.SetDestructor(&destruct_vectorlEFullTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<FullTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<FullTrack*>","std::vector<FullTrack*, std::allocator<FullTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<FullTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEFullTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<FullTrack*>*>(nullptr))->GetClass();
      vectorlEFullTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEFullTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEFullTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FullTrack*> : new vector<FullTrack*>;
   }
   static void *newArray_vectorlEFullTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FullTrack*>[nElements] : new vector<FullTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEFullTrackmUgR(void *p) {
      delete (static_cast<vector<FullTrack*>*>(p));
   }
   static void deleteArray_vectorlEFullTrackmUgR(void *p) {
      delete [] (static_cast<vector<FullTrack*>*>(p));
   }
   static void destruct_vectorlEFullTrackmUgR(void *p) {
      typedef vector<FullTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<FullTrack*>

namespace ROOT {
   static TClass *vectorlEFinalTrackVertexmUgR_Dictionary();
   static void vectorlEFinalTrackVertexmUgR_TClassManip(TClass*);
   static void *new_vectorlEFinalTrackVertexmUgR(void *p = nullptr);
   static void *newArray_vectorlEFinalTrackVertexmUgR(Long_t size, void *p);
   static void delete_vectorlEFinalTrackVertexmUgR(void *p);
   static void deleteArray_vectorlEFinalTrackVertexmUgR(void *p);
   static void destruct_vectorlEFinalTrackVertexmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<FinalTrackVertex*>*)
   {
      vector<FinalTrackVertex*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<FinalTrackVertex*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<FinalTrackVertex*>", -2, "vector", 428,
                  typeid(vector<FinalTrackVertex*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEFinalTrackVertexmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<FinalTrackVertex*>) );
      instance.SetNew(&new_vectorlEFinalTrackVertexmUgR);
      instance.SetNewArray(&newArray_vectorlEFinalTrackVertexmUgR);
      instance.SetDelete(&delete_vectorlEFinalTrackVertexmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEFinalTrackVertexmUgR);
      instance.SetDestructor(&destruct_vectorlEFinalTrackVertexmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<FinalTrackVertex*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<FinalTrackVertex*>","std::vector<FinalTrackVertex*, std::allocator<FinalTrackVertex*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<FinalTrackVertex*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEFinalTrackVertexmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<FinalTrackVertex*>*>(nullptr))->GetClass();
      vectorlEFinalTrackVertexmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEFinalTrackVertexmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEFinalTrackVertexmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FinalTrackVertex*> : new vector<FinalTrackVertex*>;
   }
   static void *newArray_vectorlEFinalTrackVertexmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FinalTrackVertex*>[nElements] : new vector<FinalTrackVertex*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEFinalTrackVertexmUgR(void *p) {
      delete (static_cast<vector<FinalTrackVertex*>*>(p));
   }
   static void deleteArray_vectorlEFinalTrackVertexmUgR(void *p) {
      delete [] (static_cast<vector<FinalTrackVertex*>*>(p));
   }
   static void destruct_vectorlEFinalTrackVertexmUgR(void *p) {
      typedef vector<FinalTrackVertex*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<FinalTrackVertex*>

namespace ROOT {
   static TClass *vectorlEFinalTrackmUgR_Dictionary();
   static void vectorlEFinalTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlEFinalTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlEFinalTrackmUgR(Long_t size, void *p);
   static void delete_vectorlEFinalTrackmUgR(void *p);
   static void deleteArray_vectorlEFinalTrackmUgR(void *p);
   static void destruct_vectorlEFinalTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<FinalTrack*>*)
   {
      vector<FinalTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<FinalTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<FinalTrack*>", -2, "vector", 428,
                  typeid(vector<FinalTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEFinalTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<FinalTrack*>) );
      instance.SetNew(&new_vectorlEFinalTrackmUgR);
      instance.SetNewArray(&newArray_vectorlEFinalTrackmUgR);
      instance.SetDelete(&delete_vectorlEFinalTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEFinalTrackmUgR);
      instance.SetDestructor(&destruct_vectorlEFinalTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<FinalTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<FinalTrack*>","std::vector<FinalTrack*, std::allocator<FinalTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<FinalTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEFinalTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<FinalTrack*>*>(nullptr))->GetClass();
      vectorlEFinalTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEFinalTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEFinalTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FinalTrack*> : new vector<FinalTrack*>;
   }
   static void *newArray_vectorlEFinalTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<FinalTrack*>[nElements] : new vector<FinalTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEFinalTrackmUgR(void *p) {
      delete (static_cast<vector<FinalTrack*>*>(p));
   }
   static void deleteArray_vectorlEFinalTrackmUgR(void *p) {
      delete [] (static_cast<vector<FinalTrack*>*>(p));
   }
   static void destruct_vectorlEFinalTrackmUgR(void *p) {
      typedef vector<FinalTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<FinalTrack*>

namespace {
  void TriggerDictionaryInitialization_InModuleTrackDict_Impl() {
    static const char* headers[] = {
"../InModuleTrack.h",
"../InModuleTrackv1.h",
"../InModuleTrackContainer.h",
"../InModuleTrackContainerv1.h",
"../FullTrack.h",
"../FullTrackv1.h",
"../FullTrackContainer.h",
"../FullTrackContainerv1.h",
"../FullTrackVertex.h",
"../FullTrackVertexv1.h",
"../FullTrackVertexContainer.h",
"../FullTrackVertexContainerv1.h",
"../FinalTrack.h",
"../FinalTrackv1.h",
"../FinalTrackContainer.h",
"../FinalTrackContainerv1.h",
"../FinalTrackVertex.h",
"../FinalTrackVertexv1.h",
"../FinalTrackVertexContainer.h",
"../FinalTrackVertexContainerv1.h",
"../TpcPolyTrack.h",
"../TpcPolyTrackv1.h",
"../TpcPolyTrackContainer.h",
"../TpcPolyTrackContainerv1.h",
"../TpcPolyCluster.h",
"../TpcPolyClusterv1.h",
"../TpcPolyClusterContainer.h",
"../TpcPolyClusterContainerv1.h",
"../TpcPolyClusterTrack.h",
"../TpcPolyClusterTrackv1.h",
"../TpcPolyClusterTrackContainer.h",
"../TpcPolyClusterTrackContainerv1.h",
"../Fitter.h",
"../IdealPadMap.h",
"../TpcPadMap.h",
"../TpcPadMapv1.h",
nullptr
    };
    static const char* includePaths[] = {
"..",
"/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/include/",
"/gpfs/mnt/gpfs02/sphenix/user/mitrankov/patterns/TPC_pattern_reco/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "InModuleTrackDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename _T1, typename _T2> struct __attribute__((annotate("$clingAutoload$bits/stl_pair.h")))  __attribute__((annotate("$clingAutoload$string")))  pair;
}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$../InModuleTrack.h")))  InModuleTrack;
class __attribute__((annotate("$clingAutoload$../InModuleTrackv1.h")))  InModuleTrackv1;
class __attribute__((annotate("$clingAutoload$../InModuleTrackContainer.h")))  InModuleTrackContainer;
class __attribute__((annotate("$clingAutoload$../InModuleTrackContainerv1.h")))  InModuleTrackContainerv1;
class __attribute__((annotate("$clingAutoload$../FullTrack.h")))  FullTrack;
class __attribute__((annotate("$clingAutoload$../FullTrackv1.h")))  FullTrackv1;
class __attribute__((annotate("$clingAutoload$../FullTrackContainer.h")))  FullTrackContainer;
class __attribute__((annotate("$clingAutoload$../FullTrackContainerv1.h")))  FullTrackContainerv1;
class __attribute__((annotate("$clingAutoload$../FullTrackVertex.h")))  FullTrackVertex;
class __attribute__((annotate("$clingAutoload$../FullTrackVertexv1.h")))  FullTrackVertexv1;
class __attribute__((annotate("$clingAutoload$../FullTrackVertexContainer.h")))  FullTrackVertexContainer;
class __attribute__((annotate("$clingAutoload$../FullTrackVertexContainerv1.h")))  FullTrackVertexContainerv1;
class __attribute__((annotate("$clingAutoload$../FinalTrack.h")))  FinalTrack;
class __attribute__((annotate("$clingAutoload$../FinalTrackv1.h")))  FinalTrackv1;
class __attribute__((annotate("$clingAutoload$../FinalTrackContainer.h")))  FinalTrackContainer;
class __attribute__((annotate("$clingAutoload$../FinalTrackContainerv1.h")))  FinalTrackContainerv1;
class __attribute__((annotate("$clingAutoload$../FinalTrackVertex.h")))  FinalTrackVertex;
class __attribute__((annotate("$clingAutoload$../FinalTrackVertexv1.h")))  FinalTrackVertexv1;
class __attribute__((annotate("$clingAutoload$../FinalTrackVertexContainer.h")))  FinalTrackVertexContainer;
class __attribute__((annotate("$clingAutoload$../FinalTrackVertexContainerv1.h")))  FinalTrackVertexContainerv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyTrack.h")))  TpcPolyTrack;
class __attribute__((annotate("$clingAutoload$../TpcPolyTrackv1.h")))  TpcPolyTrackv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyTrackContainer.h")))  TpcPolyTrackContainer;
class __attribute__((annotate("$clingAutoload$../TpcPolyTrackContainerv1.h")))  TpcPolyTrackContainerv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyCluster.h")))  TpcPolyCluster;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterv1.h")))  TpcPolyClusterv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterContainer.h")))  TpcPolyClusterContainer;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterContainerv1.h")))  TpcPolyClusterContainerv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterTrack.h")))  TpcPolyClusterTrack;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterTrackv1.h")))  TpcPolyClusterTrackv1;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterTrackContainer.h")))  TpcPolyClusterTrackContainer;
class __attribute__((annotate("$clingAutoload$../TpcPolyClusterTrackContainerv1.h")))  TpcPolyClusterTrackContainerv1;
class __attribute__((annotate("$clingAutoload$../IdealPadMap.h")))  IdealPadMap;
class __attribute__((annotate("$clingAutoload$../TpcPadMap.h")))  TpcPadMap;
class __attribute__((annotate("$clingAutoload$../TpcPadMapv1.h")))  TpcPadMapv1;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "InModuleTrackDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "../InModuleTrack.h"
#include "../InModuleTrackv1.h"
#include "../InModuleTrackContainer.h"
#include "../InModuleTrackContainerv1.h"
#include "../FullTrack.h"
#include "../FullTrackv1.h"
#include "../FullTrackContainer.h"
#include "../FullTrackContainerv1.h"
#include "../FullTrackVertex.h"
#include "../FullTrackVertexv1.h"
#include "../FullTrackVertexContainer.h"
#include "../FullTrackVertexContainerv1.h"
#include "../FinalTrack.h"
#include "../FinalTrackv1.h"
#include "../FinalTrackContainer.h"
#include "../FinalTrackContainerv1.h"
#include "../FinalTrackVertex.h"
#include "../FinalTrackVertexv1.h"
#include "../FinalTrackVertexContainer.h"
#include "../FinalTrackVertexContainerv1.h"
#include "../TpcPolyTrack.h"
#include "../TpcPolyTrackv1.h"
#include "../TpcPolyTrackContainer.h"
#include "../TpcPolyTrackContainerv1.h"
#include "../TpcPolyCluster.h"
#include "../TpcPolyClusterv1.h"
#include "../TpcPolyClusterContainer.h"
#include "../TpcPolyClusterContainerv1.h"
#include "../TpcPolyClusterTrack.h"
#include "../TpcPolyClusterTrackv1.h"
#include "../TpcPolyClusterTrackContainer.h"
#include "../TpcPolyClusterTrackContainerv1.h"
#include "../Fitter.h"
#include "../IdealPadMap.h"
#include "../TpcPadMap.h"
#include "../TpcPadMapv1.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"FinalTrack", payloadCode, "@",
"FinalTrackContainer", payloadCode, "@",
"FinalTrackContainerv1", payloadCode, "@",
"FinalTrackVertex", payloadCode, "@",
"FinalTrackVertexContainer", payloadCode, "@",
"FinalTrackVertexContainerv1", payloadCode, "@",
"FinalTrackVertexv1", payloadCode, "@",
"FinalTrackv1", payloadCode, "@",
"FullTrack", payloadCode, "@",
"FullTrackContainer", payloadCode, "@",
"FullTrackContainerv1", payloadCode, "@",
"FullTrackVertex", payloadCode, "@",
"FullTrackVertexContainer", payloadCode, "@",
"FullTrackVertexContainerv1", payloadCode, "@",
"FullTrackVertexv1", payloadCode, "@",
"FullTrackv1", payloadCode, "@",
"IdealPadMap", payloadCode, "@",
"InModuleTrack", payloadCode, "@",
"InModuleTrackContainer", payloadCode, "@",
"InModuleTrackContainerv1", payloadCode, "@",
"InModuleTrackv1", payloadCode, "@",
"TpcPadMap", payloadCode, "@",
"TpcPadMapv1", payloadCode, "@",
"TpcPolyCluster", payloadCode, "@",
"TpcPolyClusterContainer", payloadCode, "@",
"TpcPolyClusterContainerv1", payloadCode, "@",
"TpcPolyClusterTrack", payloadCode, "@",
"TpcPolyClusterTrackContainer", payloadCode, "@",
"TpcPolyClusterTrackContainerv1", payloadCode, "@",
"TpcPolyClusterTrackv1", payloadCode, "@",
"TpcPolyClusterv1", payloadCode, "@",
"TpcPolyTrack", payloadCode, "@",
"TpcPolyTrackContainer", payloadCode, "@",
"TpcPolyTrackContainerv1", payloadCode, "@",
"TpcPolyTrackv1", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("InModuleTrackDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_InModuleTrackDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_InModuleTrackDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_InModuleTrackDict() {
  TriggerDictionaryInitialization_InModuleTrackDict_Impl();
}
