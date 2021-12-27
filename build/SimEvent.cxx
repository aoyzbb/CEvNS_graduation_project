// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME SimEvent
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_SimDeposit(void *p = 0);
   static void *newArray_SimDeposit(Long_t size, void *p);
   static void delete_SimDeposit(void *p);
   static void deleteArray_SimDeposit(void *p);
   static void destruct_SimDeposit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimDeposit*)
   {
      ::SimDeposit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimDeposit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimDeposit", ::SimDeposit::Class_Version(), "SimDeposit.h", 15,
                  typeid(::SimDeposit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SimDeposit::Dictionary, isa_proxy, 4,
                  sizeof(::SimDeposit) );
      instance.SetNew(&new_SimDeposit);
      instance.SetNewArray(&newArray_SimDeposit);
      instance.SetDelete(&delete_SimDeposit);
      instance.SetDeleteArray(&deleteArray_SimDeposit);
      instance.SetDestructor(&destruct_SimDeposit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimDeposit*)
   {
      return GenerateInitInstanceLocal((::SimDeposit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimDeposit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SimTrack(void *p = 0);
   static void *newArray_SimTrack(Long_t size, void *p);
   static void delete_SimTrack(void *p);
   static void deleteArray_SimTrack(void *p);
   static void destruct_SimTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimTrack*)
   {
      ::SimTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimTrack", ::SimTrack::Class_Version(), "SimTrack.h", 16,
                  typeid(::SimTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SimTrack::Dictionary, isa_proxy, 4,
                  sizeof(::SimTrack) );
      instance.SetNew(&new_SimTrack);
      instance.SetNewArray(&newArray_SimTrack);
      instance.SetDelete(&delete_SimTrack);
      instance.SetDeleteArray(&deleteArray_SimTrack);
      instance.SetDestructor(&destruct_SimTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimTrack*)
   {
      return GenerateInitInstanceLocal((::SimTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SimEvent(void *p = 0);
   static void *newArray_SimEvent(Long_t size, void *p);
   static void delete_SimEvent(void *p);
   static void deleteArray_SimEvent(void *p);
   static void destruct_SimEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimEvent*)
   {
      ::SimEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimEvent", ::SimEvent::Class_Version(), "SimEvent.h", 17,
                  typeid(::SimEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SimEvent::Dictionary, isa_proxy, 4,
                  sizeof(::SimEvent) );
      instance.SetNew(&new_SimEvent);
      instance.SetNewArray(&newArray_SimEvent);
      instance.SetDelete(&delete_SimEvent);
      instance.SetDeleteArray(&deleteArray_SimEvent);
      instance.SetDestructor(&destruct_SimEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimEvent*)
   {
      return GenerateInitInstanceLocal((::SimEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr SimDeposit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimDeposit::Class_Name()
{
   return "SimDeposit";
}

//______________________________________________________________________________
const char *SimDeposit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimDeposit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimDeposit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimDeposit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimDeposit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimDeposit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimDeposit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimDeposit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SimTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimTrack::Class_Name()
{
   return "SimTrack";
}

//______________________________________________________________________________
const char *SimTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SimEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimEvent::Class_Name()
{
   return "SimEvent";
}

//______________________________________________________________________________
const char *SimEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void SimDeposit::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimDeposit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimDeposit::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimDeposit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimDeposit(void *p) {
      return  p ? new(p) ::SimDeposit : new ::SimDeposit;
   }
   static void *newArray_SimDeposit(Long_t nElements, void *p) {
      return p ? new(p) ::SimDeposit[nElements] : new ::SimDeposit[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimDeposit(void *p) {
      delete ((::SimDeposit*)p);
   }
   static void deleteArray_SimDeposit(void *p) {
      delete [] ((::SimDeposit*)p);
   }
   static void destruct_SimDeposit(void *p) {
      typedef ::SimDeposit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimDeposit

//______________________________________________________________________________
void SimTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimTrack(void *p) {
      return  p ? new(p) ::SimTrack : new ::SimTrack;
   }
   static void *newArray_SimTrack(Long_t nElements, void *p) {
      return p ? new(p) ::SimTrack[nElements] : new ::SimTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimTrack(void *p) {
      delete ((::SimTrack*)p);
   }
   static void deleteArray_SimTrack(void *p) {
      delete [] ((::SimTrack*)p);
   }
   static void destruct_SimTrack(void *p) {
      typedef ::SimTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimTrack

//______________________________________________________________________________
void SimEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimEvent(void *p) {
      return  p ? new(p) ::SimEvent : new ::SimEvent;
   }
   static void *newArray_SimEvent(Long_t nElements, void *p) {
      return p ? new(p) ::SimEvent[nElements] : new ::SimEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimEvent(void *p) {
      delete ((::SimEvent*)p);
   }
   static void deleteArray_SimEvent(void *p) {
      delete [] ((::SimEvent*)p);
   }
   static void destruct_SimEvent(void *p) {
      typedef ::SimEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimEvent

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 386,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace {
  void TriggerDictionaryInitialization_SimEvent_Impl() {
    static const char* headers[] = {
"/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h",
0
    };
    static const char* includePaths[] = {
"/home/mart/Software/Geant4_10.3.3/install/include/Geant4",
"/home/mart/Software/Geant4_10.3.3/install/include/Geant4",
"/home/mart/Software/ROOT_v6.22.08/install/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/G4Actions/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/PhysicsList/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/DetectorConstruction/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/ParticleGunGenerator/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/ParticleSource/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/GPSModule/include",
"/home/mart/StarXP/CDEX/cevns-ucas/source/Generator/include",
"/home/mart/Software/ROOT_v6.22.08/install/include/",
"/home/mart/StarXP/CDEX/cevns-ucas/build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "SimEvent dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$SimDeposit.h")))  __attribute__((annotate("$clingAutoload$/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h")))  SimDeposit;
class __attribute__((annotate("$clingAutoload$SimTrack.h")))  __attribute__((annotate("$clingAutoload$/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h")))  SimTrack;
class __attribute__((annotate("$clingAutoload$/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h")))  SimEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "SimEvent dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/home/mart/StarXP/CDEX/cevns-ucas/source/AnalysisManager/include/SimEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"SimDeposit", payloadCode, "@",
"SimEvent", payloadCode, "@",
"SimTrack", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("SimEvent",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_SimEvent_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_SimEvent_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_SimEvent() {
  TriggerDictionaryInitialization_SimEvent_Impl();
}
