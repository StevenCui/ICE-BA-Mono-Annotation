/******************************************************************************
 * Copyright 2017-2018 Baidu Robotic Vision Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/
#include "stdafx.h"
//#ifndef CFG_DEBUG
//#define CFG_DEBUG
//#endif
//#ifndef CFG_STEREO
//#define CFG_STEREO
//#endif
#include "IBA.h"
#include "IBA_internal.h"

//#define IBA_VERBOSE
#ifdef CFG_DEBUG
#ifdef CFG_STEREO
//#if 1
//#define IBA_DEBUG_REMOVE_RIGHT_MEASUREMENTS
#endif
//#define IBA_DEBUG_INVALIDATE_INITIAL_CAMERAS
//#define IBA_DEBUG_INVALIDATE_INITIAL_POINTS
#endif

#ifndef _WIN32
#include <boost/filesystem.hpp>
#endif

namespace IBA {

static inline void ConvertCamera(const CameraIMUState &Ci, Camera *Co) {
  if (Ci.C.R[0][0] == FLT_MAX) {
    Co->Invalidate();
    return;
  }
  Co->m_T.LA::AlignedMatrix3x3f::Set(Ci.C.R);
  Co->m_p.Set(Ci.C.p);
  Co->m_T.SetPosition(Co->m_p);
  Co->m_v.Set(Ci.v);
  Co->m_ba.Set(Ci.ba);
  Co->m_bw.Set(Ci.bw);
}

static inline void ConvertCamera(const CameraPose &Ci, Rigid3D *Co) {
  if (Ci.R[0][0] == FLT_MAX) {
    Co->Invalidate();
    return;
  }
  Co->LA::AlignedMatrix3x3f::Set(Ci.R);
  const ::Point3D p = Ci.p;
  Co->SetPosition(p);
}

static inline void ConvertCamera(const Camera &Ci, CameraIMUState *Co) {
  Ci.m_T.LA::AlignedMatrix3x3f::Get(Co->C.R);
  Ci.m_p.Get(Co->C.p);
  Ci.m_v.Get(Co->v);
  Ci.m_ba.Get(Co->ba);
  Ci.m_bw.Get(Co->bw);
}

static inline void ConvertDepth(const Depth &di, ::Depth::InverseGaussian *_do) {
  _do->u() = di.d;
  _do->s2() = di.s2;
}

static inline void ConvertDepth(const ::Depth::InverseGaussian &di, Depth *_do) {
  _do->d = di.u();
  _do->s2 = di.s2();
}

Solver::Solver() {
  m_internal = NULL;
}

Solver::~Solver() {
}

void Solver::Create(const Calibration &K, const int serial, const int verbose, const int debug,
                    const int history, const std::string param, const std::string dir,
                    const std::vector<CameraIMUState> *XsGT, const std::vector<Depth> *dsGT) {

  m_internal = new Internal();
  Camera::Calibration &_K = m_internal->m_K;
  _K.m_K.Set(K.w, K.h, K.K.fx, K.K.fy, K.K.cx, K.K.cy, K.K.ds, K.fishEye);
  const Rigid3D Tu(K.Tu);
  _K.m_Ru = Tu;
  _K.m_pu = Tu.GetTranslation();
  _K.m_ba.Set(K.ba);
  _K.m_bw.Set(K.bw);
  //_K.m_sa.Set(K.sa);
  LoadParameters(param);
  AlignedVector<Camera> &CsGT = m_internal->m_CsGT;
  if (XsGT) {
    const int N = static_cast<int>(XsGT->size());
    CsGT.Resize(N);
    for (int i = 0; i < N; ++i) {
      Camera &C = CsGT[i];
      ConvertCamera(XsGT->at(i), &C);
      C.m_ba -= _K.m_ba;
      C.m_bw -= _K.m_bw;
    }
  } else {
    CsGT.Clear();
  }
  std::vector<::Depth::InverseGaussian> &DsGT = m_internal->m_DsGT;
  if (dsGT) {
    const int N = static_cast<int>(dsGT->size());
    DsGT.resize(N);
    for (int i = 0; i < N; ++i) {
      ConvertDepth(dsGT->at(i), &DsGT[i]);
    }
  } else {
    DsGT.clear();
  }
  m_internal->m_dir = dir;
  m_internal->m_LBA.Initialize(this, serial & 255, verbose & 255,
                               static_cast<int>(static_cast<char>(debug & 255)),
                               static_cast<int>(static_cast<char>(history & 255)));
  m_internal->m_GBA.Initialize(this, (serial >> 8) & 255, (verbose >> 8) & 255,
                               static_cast<int>(static_cast<char>((debug >> 8) & 255)),
                               static_cast<int>(static_cast<char>((history >> 8) & 255)));
  m_internal->m_debug = static_cast<int>(static_cast<char>(debug & 255));
}

void Solver::Destroy() {
  delete m_internal;
  m_internal = NULL;
}

void Solver::Start() {
  m_internal->m_LBA.Start();
  m_internal->m_GBA.Start();
  m_internal->m_nFrms = 0;
  m_internal->m_Ts.resize(0);
  m_internal->m_LM.IBA_Reset();
  m_internal->m_iKF2d.assign(1, 0);
  m_internal->m_id2X.resize(0);
  m_internal->m_iX2d.resize(0);
  m_internal->m_id2idx.resize(0);
  m_internal->m_idx2iX.resize(0);
  m_internal->m_xs.resize(0);
  m_internal->m_CsLF.resize(0);
  m_internal->m_CsKF.resize(0);
  m_internal->m_ds.resize(0);
  m_internal->m_uds.resize(0);
}

void Solver::Stop() {
  m_internal->m_LBA.Stop();
  m_internal->m_GBA.Stop();
}

void Solver::PushCurrentFrame(const CurrentFrame &CF, const KeyFrame *KF, const bool serial) {
  if (KF) {
    //UT::Print("[%d]\n", KF->iFrm);
    //UT::Print("%f\n", KF->C.R[0][0]);
    if (KF->iFrm == CF.iFrm) {
      m_internal->PushCurrentFrame(CF);
      m_internal->PushKeyFrame(*KF, &m_internal->m_ILF.m_C);
      m_internal->m_LM.IBA_PushLocalFrame(LocalMap::CameraLF(m_internal->m_ILF.m_C, CF.iFrm));
      m_internal->m_LM.IBA_PushKeyFrame(m_internal->m_IKF);
      m_internal->m_LBA.PushLocalFrame(m_internal->m_ILF);
      m_internal->m_LBA.PushKeyFrame(m_internal->m_IKF, serial);
    } else {
      m_internal->PushKeyFrame(*KF);
      m_internal->PushCurrentFrame(CF);
      m_internal->m_LM.IBA_PushKeyFrame(m_internal->m_IKF);
      m_internal->m_LM.IBA_PushLocalFrame(LocalMap::CameraLF(m_internal->m_ILF.m_C, CF.iFrm));
      m_internal->m_LBA.PushKeyFrame(m_internal->m_IKF, serial);
      m_internal->m_LBA.PushLocalFrame(m_internal->m_ILF);
    }
  } else {
    m_internal->PushCurrentFrame(CF);
    m_internal->m_LM.IBA_PushLocalFrame(LocalMap::CameraLF(m_internal->m_ILF.m_C, CF.iFrm));
    m_internal->m_LBA.PushLocalFrame(m_internal->m_ILF);
  }
  m_internal->m_LBA.WakeUp(serial);
  if (m_internal->m_debug) {
    m_internal->AssertConsistency();
  }
}

void Solver::PushKeyFrame(const KeyFrame &KF, const bool serial) {
  m_internal->PushKeyFrame(KF);
  m_internal->m_LM.IBA_PushKeyFrame(m_internal->m_IKF);
  m_internal->m_LBA.PushKeyFrame(m_internal->m_IKF, serial);
  m_internal->m_LBA.WakeUp(serial);
  if (m_internal->m_debug) {
    m_internal->AssertConsistency();
  }
}

bool Solver::PushRelativeConstraint(const RelativeConstraint &Z) {
  if (Z.iFrm1 == -1 || Z.iFrm2 == -1) {
    return false;
  }
  const std::vector<LocalMap::CameraKF> &CsKF = m_internal->m_CsKF;
  const std::vector<LocalMap::CameraKF>::const_iterator i1 = std::lower_bound(CsKF.begin(),
                                                                              CsKF.end(), Z.iFrm1);
  if (i1 == CsKF.end() || i1->m_iFrm != Z.iFrm1) {
    return false;
  }
  const std::vector<LocalMap::CameraKF>::const_iterator i2 = std::lower_bound(i1, CsKF.end(),
                                                                              Z.iFrm2);
  if (i2 == CsKF.end() || i2->m_iFrm != Z.iFrm2) {
    return false;
  }
  Rigid3D T;
  LA::AlignedVector3f p;
  T.Rotation3D::Set(Z.T.R);
  p.Set(Z.T.p);
  T.SetPosition(p);
  LA::AlignedMatrix6x6f S;
  S.Set(Z.S.S);
  CameraPrior::Pose _Z;
  if (!_Z.Initialize(BA_WEIGHT_PRIOR_CAMERA_RELATIVE_CONSTRAINT,
                     static_cast<int>(i1 - CsKF.begin()),
                     static_cast<int>(i2 - CsKF.begin()), T, S)) {
    return false;
  }
  //////////////////////////////////////////////////////////////////////////
  //const int verboseBkp = m_internal->m_GBA.m_verbose;
  //m_internal->m_GBA.m_verbose = 2;
  //////////////////////////////////////////////////////////////////////////
  m_internal->m_GBA.PushCameraPriorPose(Z.iFrm2, _Z);
  m_internal->m_GBA.WakeUp();
  //////////////////////////////////////////////////////////////////////////
  //m_internal->m_GBA.m_verbose = verboseBkp;
  //////////////////////////////////////////////////////////////////////////
  return true;
}

void Solver::SetCallbackLBA(const IbaCallback& iba_callback) {
  m_internal->m_LBA.SetCallback(iba_callback);
}

void Solver::SetCallbackGBA(const IbaCallback& iba_callback) {
  m_internal->m_GBA.SetCallback(iba_callback);
}

bool Solver::GetSlidingWindow(SlidingWindow *SW) {
  const ubyte Uc = m_internal->m_LM.IBA_Synchronize(m_internal->m_nFrms,
                                                    m_internal->m_CsLF, m_internal->m_CsKF,
                                                    m_internal->m_ds, m_internal->m_uds);
  SW->iFrms.resize(0);
  SW->CsLF.resize(0);
  if (Uc & LM_FLAG_FRAME_UPDATE_CAMERA_LF) {
    CameraIMUState C;
    for (std::list<LocalMap::CameraLF>::const_iterator i = m_internal->m_CsLF.begin();
         i != m_internal->m_CsLF.end(); ++i) {
      if (!i->m_uc) {
        continue;
      }
      SW->iFrms.push_back(i->m_iFrm);
      ConvertCamera(i->m_C, &C);
      SW->CsLF.push_back(C);
    }
  }
  SW->iFrmsKF.resize(0);
  SW->CsKF.resize(0);
  const std::vector<LocalMap::CameraKF> &CsKF = m_internal->m_CsKF;
  if (Uc & LM_FLAG_FRAME_UPDATE_CAMERA_KF) {
    ::Point3D p;
    Rigid3D C1;
    CameraPose C2;
    const int nKFs = static_cast<int>(CsKF.size());
    for (int iKF = 0; iKF < nKFs; ++iKF) {
      const LocalMap::CameraKF &C = CsKF[iKF];
      if (!(C.m_uc & LM_FLAG_FRAME_UPDATE_CAMERA_KF)) {
        continue;
      }
      C1 = C.m_C;
      C1.LA::AlignedMatrix3x3f::Get(C2.R);
      C1.GetPosition(p);
      p.Get(C2.p);
      SW->iFrmsKF.push_back(C.m_iFrm);
      SW->CsKF.push_back(C2);
    }
  }
  SW->Xs.resize(0);
  const ubyte ucFlag = LM_FLAG_FRAME_UPDATE_CAMERA_KF | LM_FLAG_FRAME_UPDATE_DEPTH;
  //if (Uc & LM_FLAG_FRAME_UPDATE_DEPTH) {
  if (Uc & ucFlag) {
    Rigid3D C;
    ::Point3D X1;
    Point3D X2;
    const int nKFs = static_cast<int>(CsKF.size());
    for (int iKF = 0; iKF < nKFs; ++iKF) {
      const LocalMap::CameraKF &_C = CsKF[iKF];
      //if (!(_C.m_uc & LM_FLAG_FRAME_UPDATE_DEPTH)) {
      if (!(_C.m_uc & ucFlag)) {
        continue;
      }
      const bool uc = (_C.m_uc & LM_FLAG_FRAME_UPDATE_CAMERA_KF) != 0;
      const int id1 = m_internal->m_iKF2d[iKF], id2 = m_internal->m_iKF2d[iKF + 1];
      const ::Point2D *xs = m_internal->m_xs.data() + id1;
      const ::Depth::InverseGaussian *ds = m_internal->m_ds.data() + id1;
      const ubyte *uds = m_internal->m_uds.data() + id1;
      const int *idxs = m_internal->m_id2idx.data() + id1;
      const int Nx = id2 - id1;
      C = _C.m_C;
      for (int ix = 0; ix < Nx; ++ix) {
        //if (!(uds[ix] & LM_FLAG_TRACK_UPDATE_DEPTH)) {
        if (!uc && !(uds[ix] & LM_FLAG_TRACK_UPDATE_DEPTH)) {
          continue;
        }
        X2.idx = idxs[ix];
        if (X2.idx == -1) {
          continue;
        }
        C.ApplyInversely(xs[ix], 1.0f / ds[ix].u(), X1);
        X1.Get(X2.X);
        SW->Xs.push_back(X2);
      }
    }
  }
  return Uc != LM_FLAG_FRAME_DEFAULT;
}

void Solver::PrintSlidingWindow(const SlidingWindow &SW) {
  UT::Print("\n");
  const int nLFs = static_cast<int>(SW.iFrms.size());
  UT::Print("LF = %d {", nLFs);
  for (int i = 0; i < nLFs; ++i) {
    UT::Print(" %d", SW.iFrms[i]);
  }
  UT::Print("}\n");
  const int nKFs = static_cast<int>(SW.iFrmsKF.size());
  UT::Print("KF = %d {", nKFs);
  for (int i = 0; i < nKFs; ++i) {
    UT::Print(" %d", SW.iFrmsKF[i]);
  }
  UT::Print("}\n");
  const int NX = static_cast<int>(SW.Xs.size());
  UT::Print(" X = %d {", NX);
  for (int i = 0; i < NX; ++i) {
    UT::Print(" %d", SW.Xs[i].idx);
  }
  UT::Print("}\n");
}

int Solver::GetKeyFrames() {
  return static_cast<int>(m_internal->m_CsKF.size());
}

int Solver::GetKeyFrameIndex(const int iKF) {
  return m_internal->m_CsKF[iKF].m_iFrm;
}

int Solver::SearchKeyFrame(const int iFrm) {
  const std::vector<LocalMap::CameraKF> &CsKF = m_internal->m_CsKF;
  const std::vector<LocalMap::CameraKF>::const_iterator i = std::lower_bound(CsKF.begin(),
                                                                             CsKF.end(), iFrm);
  if (i == CsKF.end() || i->m_iFrm != iFrm) {
    return -1;
  }
  return static_cast<int>(i - CsKF.begin());
}

bool Solver::DeleteKeyFrame(const int iFrm) {
  const int iKF = SearchKeyFrame(iFrm);
  if (iKF == -1) {
    return false;
  }
  std::vector<LocalMap::CameraKF> &CsKF = m_internal->m_CsKF;
  const int nKFs = static_cast<int>(CsKF.size());
  if (iKF == nKFs - 1) {
    return false;
  }

  CsKF.erase(CsKF.begin() + iKF);

  std::vector<int> &iKF2d = m_internal->m_iKF2d;
  const int id1 = iKF2d[iKF], id2 = iKF2d[iKF + 1], Nd = id2 - id1;
  for (int jKF = iKF + 1; jKF <= nKFs; ++jKF) {
    iKF2d[jKF] -= Nd;
  }
  iKF2d.erase(iKF2d.begin() + iKF);
  std::vector<int> &id2X = m_internal->m_id2X, &iX2d = m_internal->m_iX2d;
  std::vector<int> &id2idx = m_internal->m_id2idx, &idx2iX = m_internal->m_idx2iX;
  for (int id = id1; id < id2; ++id) {
    const int iX = id2X[id], idx = id2idx[id];
    if (iX != -1) {
      iX2d[iX] = idx2iX[idx] = -1;
    }
  }
  const int ND = static_cast<int>(id2X.size());
  for (int id = id2; id < ND; ++id) {
    const int iX = id2X[id];
    if (iX != -1) {
      iX2d[iX] -= Nd;
    }
  }
  id2X.erase(id2X.begin() + id1, id2X.begin() + id2);
  id2idx.erase(id2idx.begin() + id1, id2idx.begin() + id2);
  std::vector<::Point2D> &xs = m_internal->m_xs;
  xs.erase(xs.begin() + id1, xs.begin() + id2);
  std::vector<::Depth::InverseGaussian> &ds = m_internal->m_ds;
  ds.erase(ds.begin() + id1, ds.begin() + id2);
  std::vector<ubyte> &uds = m_internal->m_uds;
  uds.erase(uds.begin() + id1, uds.begin() + id2);
  m_internal->m_LM.IBA_DeleteKeyFrame(iFrm, iKF);
  m_internal->m_LBA.PushDeleteKeyFrame(iFrm, iKF);
  m_internal->m_LBA.WakeUp();
  return true;
}

void Solver::GetMapPointIndexes(std::vector<int> *idxs) {
  idxs->resize(0);
  const std::vector<int> &idx2iX = m_internal->m_idx2iX;
  const int N = static_cast<int>(idx2iX.size());
  for (int idx = 0; idx < N; ++idx) {
    if (idx2iX[idx] != -1) {
      idxs->push_back(idx);
    }
  }
}

void Solver::DeleteMapPoints(const std::vector<int> &idxs) {
  if (idxs.empty()) {
    return;
  }
  std::vector<int> &ids = m_internal->m_idxsTmp;
  std::vector<int> &id2X = m_internal->m_id2X, &iX2d = m_internal->m_iX2d;
  std::vector<int> &id2idx = m_internal->m_id2idx, &idx2iX = m_internal->m_idx2iX;
  const int N = static_cast<int>(idxs.size());
  ids.resize(0);
  for (int i = 0; i < N; ++i) {
    const int idx = idxs[i], iX = idx2iX[idx];
    if (iX == -1) {
      continue;
    }
    const int id = iX2d[iX];
    ids.push_back(id);
    id2X[id] = iX2d[iX] = id2idx[id] = idx2iX[idx] = -1;
  }
  std::sort(ids.begin(), ids.end());

  m_internal->m_LBA.PushDeleteMapPoints(m_internal->m_nFrms, ids);
}

void Solver::UpdateCameras(const std::vector<int> &iFrms, const std::vector<CameraPose> &Cs,
                           const bool serial) {
  Rigid3D C;
  std::vector<GlobalMap::InputCamera> &ICs = m_internal->m_ICs;
  ICs.resize(0);
  const int N = static_cast<int>(iFrms.size());

  for (int i = 0; i < N; ++i) {
    ConvertCamera(Cs[i], &C);
    ICs.push_back(GlobalMap::InputCamera(C, iFrms[i]));
  }
  if (ICs.empty()) {
    return;
  }

  m_internal->m_LBA.PushUpdateCameras(m_internal->m_nFrms, ICs, serial);
  m_internal->m_LBA.WakeUp(serial);
}

bool Solver::PropagateState(const std::vector<IMUMeasurement> &us, const float t1, const float t2,
                            const CameraIMUState &X1, CameraIMUState *X2, CameraPoseCovariance *S) {
  AlignedVector<IMU::Measurement> &_us = m_internal->m_us;
  const int N = static_cast<int>(us.size());
  _us.Resize(N);
  const Rotation3D Ru = m_internal->m_K.m_Ru;
  const ::Point3D pu = m_internal->m_K.m_pu;
  const LA::AlignedVector3f ba = m_internal->m_K.m_ba;
  const LA::AlignedVector3f bw = m_internal->m_K.m_bw;
  for (int i = 0; i < N; ++i) {
    const IMUMeasurement &u = us[i];
    IMU::Measurement &_u = _us[i];
    _u.m_a.Set(u.a);
    _u.m_w.Set(u.w);
    if (Ru.Valid()) {
      _u.m_a = Ru.GetApplied(_u.m_a);
      _u.m_w = Ru.GetApplied(_u.m_w);
    }
    _u.m_a -= ba;
    _u.m_w -= bw;
    _u.t() = u.t;
  }
  Camera C1, C2;
  IMU::Delta &D = m_internal->m_D;
  ConvertCamera(X1, &C1);
  IMU::PreIntegrate(_us, t1, t2, C1, &D, &m_internal->m_work, S != NULL, NULL, NULL,
                    BA_ANGLE_EPSILON);
  IMU::Propagate(pu, D, C1, C2, BA_ANGLE_EPSILON);
  ConvertCamera(C2, X2);
  if (S != NULL) {
    LA::AlignedMatrixXf W;
    W.Resize(15, 15);
    for (int i = 0, ip = 0; i < 5; ++i, ip += 3) {
      for (int j = i, jp = ip; j < 5; ++j, jp += 3) {
        W.SetBlock(ip, jp, D.m_W[i][j]);
      }
    }
    if (!W.InverseLDL()) {
      return false;
    }
    LA::AlignedMatrix3x3f Srr, Srp, Spr, Spp;
    W.GetBlock(0, 0, Srr);  W.GetBlock(0, 6, Srp);
    W.GetBlock(6, 0, Spr);  W.GetBlock(6, 6, Spp);
    LA::AlignedMatrix3x3f U, Upp, Upr/*, Urp*/, Urr;
    const SkewSymmetricMatrix u = D.m_RT.GetApplied(pu);
    const Rotation3D R1T = C1.m_T.GetRotationTranspose();
    SkewSymmetricMatrix::ABT(Spr, u, U);
    U += Spp;
    LA::AlignedMatrix3x3f::ABT(R1T, U, Upp);
    SkewSymmetricMatrix::ABT(Srr, u, U);
    U += Srp;
    LA::AlignedMatrix3x3f::ABT(R1T, U, Upr);
    //LA::AlignedMatrix3x3f::ABT(R1T, Spr, Urp);
    LA::AlignedMatrix3x3f::ABT(R1T, Srr, Urr);
    SkewSymmetricMatrix::ABT(Upr, u, U);
    U += Upp;
    LA::AlignedMatrix3x3f::ABT(U, R1T, Spp);
    LA::AlignedMatrix3x3f::ABT(Upr, R1T, Spr);
    LA::AlignedMatrix3x3f::ABT(Urr, R1T, Srr);
    Spr.GetTranspose(Srp);
    LA::AlignedMatrix6x6f _S;
    _S.Set(Spp, Spr, Srp, Srr);
    _S.Get(S->S);
  }
  return true;
}

bool Solver::SaveFeatures(const std::string fileName, const std::vector<MapPointMeasurement> &zs,
                          const std::vector<MapPoint> *Xs) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  fprintf(fp, "key_points: [ ");
  const int Nz = static_cast<int>(zs.size());
  for (int iz = 0; iz < Nz; ++iz) {
    const MapPointMeasurement &z = zs[iz];
    //const int id = m_internal->m_iX2d[z.iX];
    const int iX = m_internal->m_idx2iX[z.idx];
    const int id = m_internal->m_iX2d[iX];
    if (id != -1) {
      fprintf(fp, "%f, %f, 0, 0, 0, 0, %d,\n", z.x.x[0], z.x.x[1], m_internal->m_id2idx[id]);
    }
  }
  fprintf(fp, "]\n");
  if (Xs) {
    fprintf(fp, "map_points: [\n");
    const int Nx = static_cast<int>(Xs->size());
    fprintf(fp, "%d\n", Nx);
    for (int ix = 0; ix < Nx; ++ix) {
      const MapPoint &X = Xs->at(ix);
      const int Nz = static_cast<int>(X.zs.size());
      fprintf(fp, "%d %d\n", X.X.idx, Nz);
      for (int i = 0; i < Nz; ++i) {
        const MapPointMeasurement &z = X.zs[i];
        fprintf(fp, "%d %f %f\n", z.iFrm, z.x.x[0], z.x.x[1]);
      }
    }
    fprintf(fp, "]\n");
  }
  fclose(fp);
  return true;
}

bool Solver::LoadFeatures(const std::string fileName, const int iFrm,
                          std::vector<MapPointMeasurement> *zs,
                          std::vector<MapPoint> *Xs,
                          const std::vector<int> *iFrms) {
  {
    zs->resize(0);
    if (Xs) {
      Xs->resize(0);
    }
  }
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    return false;
  }
  char line[UT_STRING_WIDTH_MAX];
  MapPointMeasurement z;
  int i = 0/*, idx*/;
  char *num;
  float v;
  const bool cov = FTR_VARIANCE == 0.0f;
  if (!cov) {
    z.x.S[0][0] = z.x.S[1][1] = FTR_VARIANCE;
    z.x.S[0][1] = z.x.S[1][0] = 0.0f;
  }
  const std::vector<int> &idx2iX = m_internal->m_idx2iX;
  const int idxMax = static_cast<int>(idx2iX.size()) - 1;
  while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
    const int len = int(std::remove(line, std::remove(line, line + strlen(line), 10), ' ') - line);
    if (len == 0 || line[0] == '%' || line[len - 1] == ':') {
      continue;
    }
    line[len] = 0;
    num = strtok(line, " ,");
    while (num && (sscanf(num, "%f", &v) == 1 || sscanf(num, "key_points:[%f", &v) == 1)) {
      switch (i) {
      case 0:
        z.x.x[0] = v;
        break;
      case 1:
        z.x.x[1] = v;
        break;
      case 2:
        if (cov) {
          z.x.S[0][0] = v;
        }
        break;
      case 3:
        if (cov) {
          z.x.S[0][1] = z.x.S[0][1] = v;
        }
        break;
      case 4:
        if (cov) {
          z.x.S[1][1] = v;
        }
        break;
      case 5:
        break;
      case 6:
        //idx = static_cast<int>(v);
        //z.iX = idx > idxMax ? -1 : idx2iX[idx];
        z.idx = static_cast<int>(v);
        break;
      }
      ++i;
      if (i == 7) {
        i = 0;
        if (z.idx > idxMax || idx2iX[z.idx] == -1) {
        //if (z.iX == -1) {
          if (Xs) {
            {
              const int iX = static_cast<int>(Xs->size());
              Xs->resize(iX + 1);
              MapPoint &X = Xs->back();
              X.X.idx = z.idx;
              X.X.X[0] = FLT_MAX;
              z.iFrm = iFrm;
              X.zs.assign(1, z);
            }
          }
        } else {
          zs->push_back(z);
        }
      }
      num = strtok(NULL, " ,");
    }
    if (line[len - 1] == ']') {
      break;
    }
  }
  if (Xs) {
    int Nx, Nz, iX, idx;
    float tmp[2];
    while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
      const int len = int(std::remove(line, std::remove(line, line + strlen(line), 10), ' ') - line);
      if (len == 0 || line[0] == '%' || line[len - 1] == ':') {
        continue;
      }
      line[len] = 0;
      if (strcmp(line, "map_points:[") != 0) {
        continue;
      }
      fscanf(fp, "%d", &Nx);
      {
        Xs->reserve(Xs->size() + Nx);
      }
      for (int ix = 0; ix < Nx; ++ix) {
        fscanf(fp, "%d", &idx);
        {
          iX = static_cast<int>(Xs->size());
          Xs->resize(iX + 1);
          Xs->back().zs.resize(0);
        }
        MapPoint &X = Xs->at(iX);
        X.X.idx = idx;
        fscanf(fp, "%d", &Nz);
        X.X.X[0] = FLT_MAX;
        for (int i = 0; i < Nz; ++i) {
          fscanf(fp, "%d %f %f", &z.iFrm, &z.x.x[0], &z.x.x[1]);
          z.x.S[0][0] = z.x.S[1][1] = FTR_VARIANCE;
          z.x.S[0][1] = z.x.S[1][0] = 0.0f;
          if (iFrms) {
            z.iFrm = iFrms->at(z.iFrm);
            if (z.iFrm == -1) {
              continue;
            }
          }
          X.zs.insert(std::lower_bound(X.zs.begin(), X.zs.end(), z), z);
        }
      }
      break;
    }
  }
  fclose(fp);
  return true;
}

void Solver::SaveB(FILE *fp) {
  m_internal->m_LM.SaveB(fp);
  m_internal->m_GM.SaveB(fp);
  m_internal->m_LBA.SaveB(fp);
  m_internal->m_GBA.SaveB(fp);
  UT::SaveB(m_internal->m_nFrms, fp);
  FRM::VectorSaveB(m_internal->m_Ts, fp);
  UT::VectorSaveB(m_internal->m_iKF2d, fp);
  UT::VectorSaveB(m_internal->m_id2X, fp);
  UT::VectorSaveB(m_internal->m_iX2d, fp);
  UT::VectorSaveB(m_internal->m_id2idx, fp);
  UT::VectorSaveB(m_internal->m_idx2iX, fp);
  UT::VectorSaveB(m_internal->m_xs, fp);
  UT::ListSaveB(m_internal->m_CsLF, fp);
  UT::VectorSaveB(m_internal->m_CsKF, fp);
  UT::VectorSaveB(m_internal->m_ds, fp);
  UT::VectorSaveB(m_internal->m_uds, fp);
}

void Solver::LoadB(FILE *fp) {
  m_internal->m_LM.LoadB(fp);
  m_internal->m_GM.LoadB(fp);
  m_internal->m_LBA.LoadB(fp);
  m_internal->m_GBA.LoadB(fp);
  UT::LoadB(m_internal->m_nFrms, fp);
  FRM::VectorLoadB(m_internal->m_Ts, fp);
  UT::VectorLoadB(m_internal->m_iKF2d, fp);
  UT::VectorLoadB(m_internal->m_id2X, fp);
  UT::VectorLoadB(m_internal->m_iX2d, fp);
  UT::VectorLoadB(m_internal->m_id2idx, fp);
  UT::VectorLoadB(m_internal->m_idx2iX, fp);
  UT::VectorLoadB(m_internal->m_xs, fp);
  UT::ListLoadB(m_internal->m_CsLF, fp);
  UT::VectorLoadB(m_internal->m_CsKF, fp);
  UT::VectorLoadB(m_internal->m_ds, fp);
  UT::VectorLoadB(m_internal->m_uds, fp);
}

void Solver::ComputeErrorLBA(Error *e) {
  m_internal->m_LBA.ComputeErrorFeature(&e->ex);
  m_internal->m_LBA.ComputeErrorIMU(&e->eur, &e->eup, &e->euv, &e->euba, &e->eubw);
  m_internal->m_LBA.ComputeErrorDrift(&e->edr, &e->edp);
}

void Solver::ComputeErrorGBA(Error *e) {
  m_internal->m_GBA.ComputeErrorFeature(&e->ex);
  m_internal->m_GBA.ComputeErrorIMU(&e->eur, &e->eup, &e->euv, &e->euba, &e->eubw);
  m_internal->m_GBA.ComputeErrorDrift(&e->edr, &e->edp);
}

float Solver::ComputeRMSELBA() {
  return m_internal->m_LBA.ComputeRMSE();
}

float Solver::ComputeRMSEGBA() {
  return m_internal->m_GBA.ComputeRMSE();
}

float Solver::GetTotalDistance() {
  return m_internal->m_LBA.GetTotalDistance();
}

static inline std::string ConvertFileName(const std::string fileName, const std::string dir, const bool append = true) {
  if (fileName == "") {
    return "";
  }
  std::string _fileName;
  if (UT::FileNameRemoveDirectoryExtension(fileName) == "*") {
    _fileName = dir.substr(0, dir.length() - 1) + "." + UT::FileNameExtractExtension(fileName);
    _fileName = UT::FileNameReplaceDirectory(_fileName, UT::FileNameExtractDirectory(_fileName),
                                             UT::FileNameExtractDirectory(fileName));
  } else {
    _fileName = fileName;
  }
  _fileName = UT::FileNameReplaceDirectory(_fileName, ".", dir);
  if (append) {
    _fileName = UT::FileNameAppendSuffix(_fileName);
  }
  const std::string _dir = UT::FileNameExtractDirectory(_fileName);
  if (!UT::FileExists(_dir) && _dir != "") {
#ifdef WIN32
    CreateDirectory(_dir.c_str(), 0);
#else
    boost::filesystem::create_directories(_dir);
#endif
  }
  return _fileName;
}

bool Solver::SaveCamerasLBA(const std::string fileName, const bool append, const bool poseOnly) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SaveCameras(_fileName, poseOnly);
}

bool Solver::SaveCamerasGBA(const std::string fileName, const bool append, const bool poseOnly) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_GBA.SaveCameras(_fileName, poseOnly);
}

bool Solver::SaveCostsLBA(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SaveCosts(_fileName, type);
}

bool Solver::SaveCostsGBA(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_GBA.SaveCosts(_fileName, type);
}

bool Solver::SaveResidualsLBA(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SaveResiduals(_fileName, type);
}

bool Solver::SaveResidualsGBA(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_GBA.SaveResiduals(_fileName, type);
}

bool Solver::SavePriors(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SavePriors(_fileName, type);
}

bool Solver::SaveMarginalizations(const std::string fileName, const bool append, const int type) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SaveMarginalizations(_fileName, type);
}

bool Solver::SavePointsLBA(const std::string fileName, const bool append) {
  return m_internal->SavePoints(m_internal->m_LBA.m_CsKF, m_internal->m_LBA.m_ds,
                                fileName, append);
}

bool Solver::SavePointsGBA(const std::string fileName, const bool append) {
  return m_internal->SavePoints(m_internal->m_GBA.m_Cs, m_internal->m_GBA.m_ds,
                                fileName, append);
}

bool Internal::SavePoints(const AlignedVector<Rigid3D> &CsKF,
                          const std::vector<::Depth::InverseGaussian> &ds,
                          const std::string fileName, const bool append) {
  const std::string _fileName = ConvertFileName(fileName, m_dir, append);
  FILE *fp = fopen(_fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  AlignedVector<::Point3D> Xs;
  const int Nd = static_cast<int>(ds.size());
  m_work.Resize(Xs.BindSize(Nd));
  Xs.Bind(m_work.Data(), Nd);
  const int nKFs = CsKF.Size();
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const Rigid3D &C = CsKF[iKF];
    const int id1 = m_iKF2d[iKF], id2 = m_iKF2d[iKF + 1];
    const ::Point2D *xs = m_xs.data() + id1;
    const ::Depth::InverseGaussian *_ds = ds.data() + id1;
    const int *idxs = m_id2idx.data() + id1;
    ::Point3D *_Xs = Xs.Data() + id1;
    const int Nx = id2 - id1;
    for (int ix = 0; ix < Nx; ++ix) {
      if (idxs[ix] != -1) {
        C.ApplyInversely(xs[ix], 1.0f / _ds[ix].u(), _Xs[ix]);
      }
    }
  }
  const int N = static_cast<int>(m_idx2iX.size());
  for (int idx = 0; idx < N; ++idx) {
    const int iX = m_idx2iX[idx];
    if (iX == -1) {
      continue;
    }
    const int id = m_iX2d[iX];
    if (id == -1) {
      continue;
    }
    const ::Point3D &X = Xs[id];
    fprintf(fp, "%d %f %f %f\n", idx, X.x(), X.y(), X.z());
  }
  fclose(fp);
  UT::PrintSaved(_fileName);
  return true;
}

void Solver::GetTimeLBA(Time *t) {
  t->t = m_internal->m_LBA.GetTotalTime(&t->n);
  //t->nd = m_internal->m_LBA.GetTotalPoints();
}

void Solver::GetTimeGBA(Time *t) {
  t->t = m_internal->m_GBA.GetTotalTime(&t->n);
  //t->nd = m_internal->m_GBA.GetTotalPoints();
}

bool Solver::SaveTimesLBA(const std::string fileName, const bool append) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_LBA.SaveTimes(_fileName);
}

bool Solver::SaveTimesGBA(const std::string fileName, const bool append) {
  const std::string _fileName = ConvertFileName(fileName, m_internal->m_dir, append);
  return m_internal->m_GBA.SaveTimes(_fileName);
}

static inline void Rectify(const ::Intrinsic &K, const Point2D &x,
                           ::Point2D *xn, LA::SymmetricMatrix2x2f *Wn,
                           ::Intrinsic::UndistortionMap *UM = NULL,
                           const float eps = 0.0f) {
  ::Point2D xd;
  Point2DCovariance Sd, Sn;
  LA::AlignedMatrix2x2f J, JS;
  LA::AlignedMatrix2x2f Wd;
  K.k().ImageToNormalized(x.x, xd);
  K.Undistort(xd, xn, &J, UM);

  K.k().ImageToNormalized(x.S, Sd);
  Sd.GetInverse(Wd);
  LA::AlignedMatrix2x2f::ABT(J, Wd, JS);
  LA::SymmetricMatrix2x2f::ABT(JS, J, *Wn);
  Wn->GetInverse(Sn);

  Sn.IncreaseDiagonal(eps);
  Sn.GetInverse(*Wn);
  //////////////////////////////////////////////////////////////////////////
  //Wn->Set(Wd.m00(), Wd.m01(), Wd.m11());
  //////////////////////////////////////////////////////////////////////////
}

const LocalBundleAdjustor::InputLocalFrame& Internal::PushCurrentFrame(const CurrentFrame &CF) {
  m_ILF.m_T.m_iFrm = CF.iFrm;
  m_ILF.m_T.m_t = CF.t;
  m_ILF.m_T.m_fileName = CF.fileName;

  m_Ts.push_back(m_ILF.m_T);
  ConvertCamera(CF.C, &m_ILF.m_C);
  ConvertDepth(CF.d, &m_ILF.m_d);
  ConvertFeatureMeasurements(CF.zs, &m_ILF);
  const int Nu = static_cast<int>(CF.us.size());
  m_ILF.m_us.Resize(Nu);
  for (int i = 0; i < Nu; ++i) {
    const IMUMeasurement &u = CF.us[i];
    IMU::Measurement &_u = m_ILF.m_us[i];
    _u.m_a.Set(u.a);
    _u.m_w.Set(u.w);
    if (m_K.m_Ru.Valid()) {
      _u.m_a = m_K.m_Ru.GetApplied(_u.m_a);
      _u.m_w = m_K.m_Ru.GetApplied(_u.m_w);
    }
    _u.m_a -= m_K.m_ba;
    _u.m_w -= m_K.m_bw;
    _u.t() = u.t;
  }
  ++m_nFrms;
  return m_ILF;
}

const GlobalMap::InputKeyFrame& Internal::PushKeyFrame(const KeyFrame &KF, const Camera *C) {

  const std::vector<FRM::Tag>::iterator T = std::lower_bound(m_Ts.begin(), m_Ts.end(), KF.iFrm);

  m_IKF.m_T = *T;
  m_Ts.erase(m_Ts.begin(), T);
  ConvertCamera(KF.C, &m_IKF.m_C.m_T);
  m_IKF.m_C.m_p.Set(KF.C.p);
  if (C) {
    m_IKF.m_C.m_v = C->m_v;
    m_IKF.m_C.m_ba = C->m_ba;
    m_IKF.m_C.m_bw = C->m_bw;
  } else {
    m_IKF.m_C.m_v.Invalidate();
    m_IKF.m_C.m_ba.Invalidate();
    m_IKF.m_C.m_bw.Invalidate();
  }
  ConvertDepth(KF.d, &m_IKF.m_d);
  ConvertFeatureMeasurements(KF.zs, &m_IKF);
  m_CsKF.push_back(LocalMap::CameraKF(m_IKF.m_C.m_T, KF.iFrm));

  const int Nx1 = static_cast<int>(KF.Xs.size());
   }
  }

  int idxMax = static_cast<int>(m_idx2iX.size()) - 1;
  std::vector<MapPointIndex> &ixsSort = m_idxsSortTmp;
  ixsSort.resize(Nx1);
  for (int ix = 0; ix < Nx1; ++ix) {
    const MapPoint &X = KF.Xs[ix];
    idxMax = std::max(idxMax, X.X.idx);
    ixsSort[ix].Set(X.zs.front().iFrm, X.X.idx, ix);
  }
  std::sort(ixsSort.begin(), ixsSort.end());
  int iKF0, iKF1, iKF2;
  Rigid3D C1;
  if (ixsSort.empty()) {
    iKF0 = iKF2 = -1;
  } else {
    LocalMap::CameraKF _C;
    if ((_C.m_iFrm = ixsSort.front().m_iFrm) < m_CsKF.front().m_iFrm) {
      iKF0 = 0;
    } else {
      iKF0 = static_cast<int>(std::upper_bound(m_CsKF.begin(), m_CsKF.end(), _C) -
                                               m_CsKF.begin()) - 1;
    }
    _C.m_iFrm = ixsSort.back().m_iFrm;
    iKF2 = static_cast<int>(std::upper_bound(m_CsKF.begin() + iKF0, m_CsKF.end(), _C) -
                                             m_CsKF.begin());
    C1 = m_CsKF[iKF0].m_C;
  }
  int i1, i2;
  FTR::Measurement z;

  ::Point3D X;
  m_IKF.m_Xs.resize(Nx1);
  std::vector<int> &ixs = m_idxsTmp;
  ixs.resize(Nx1);
  const float eps = FTR_VARIANCE_EPSILON * m_K.m_K.fxyI();
  const int nKFs = static_cast<int>(m_CsKF.size());
  for (i1 = i2 = 0, iKF1 = iKF0; i1 < Nx1; ++i1) {
    const MapPointIndex &ix = ixsSort[i1];
    if (ix.m_iFrm < m_CsKF[iKF1].m_iFrm) {
      continue;
    } else if (ix.m_iFrm > m_CsKF[iKF1].m_iFrm) {
      iKF1 = static_cast<int>(std::lower_bound(m_CsKF.begin() + iKF1,
                                               m_CsKF.begin() + iKF2, ix.m_iFrm) -
                                               m_CsKF.begin());
      if (iKF1 == nKFs || m_CsKF[iKF1].m_iFrm != ix.m_iFrm) {
        --iKF1;
        continue;
      }
      C1 = m_CsKF[iKF1].m_C;
    }
    z.m_iKF = iKF1;
    const MapPoint &X1 = KF.Xs[ix.m_ix];
    GlobalMap::Point &X2 = m_IKF.m_Xs[i2];
    X2.m_zs.resize(0);
    const int Nz1 = static_cast<int>(X1.zs.size());
    for (int i = 0; i < Nz1; ++i) {
      const MapPointMeasurement &_z = X1.zs[i];
      {
        Rectify(m_K.m_K, _z.x, &z.m_z, &z.m_W, &m_UM, eps);
        {
          z.m_iKF = static_cast<int>(std::lower_bound(m_CsKF.begin() + z.m_iKF,
                                                      m_CsKF.begin() + iKF2, _z.iFrm) -
                                                      m_CsKF.begin());
          if (z.m_iKF != nKFs && m_CsKF[z.m_iKF].m_iFrm == _z.iFrm) {
            X2.m_zs.push_back(z);
          } else {
            --z.m_iKF;
          }
        }
      }
    }
    const int Nz2 = static_cast<int>(X2.m_zs.size());
    if (Nz2 == 0) {
      continue;
    }
    const FTR::Measurement &_z = X2.m_zs.front();
    X2.m_iKF = _z.m_iKF;
    X2.m_x.m_x = _z.m_z;
    X2.m_W = _z.m_W;
    X2.m_zs.erase(X2.m_zs.begin());
    X.Set(X1.X.X);
    if (X.Valid()) {
      X2.m_d.Initialize(1.0f / C1.GetAppliedZ(X));
    } else {
      X2.m_d.Invalidate();
    }
    ixs[i2] = ix.m_ix;
    ++i2;
  }
  const int Nx2 = i2;
  m_IKF.m_Xs.resize(Nx2);
  ixs.resize(Nx2);

  m_iKF2d.push_back(m_iKF2d.back());
  const int NX1 = static_cast<int>(m_iX2d.size()), NX2 = NX1 + Nx1;
  m_iX2d.resize(NX2, -1);
  m_idx2iX.resize(idxMax + 1, -1);
  for (int i1 = 0, i2 = 0; i1 < Nx2; i1 = i2) {
    const int iKF = m_IKF.m_Xs[i1].m_iKF;
    for (i2 = i1 + 1; i2 < Nx2 && m_IKF.m_Xs[i2].m_iKF == iKF; ++i2) {}
    const int id = m_iKF2d[iKF + 1], Nx = i2 - i1;
    for (int jKF = iKF; jKF < nKFs; ++jKF) {
      m_iKF2d[jKF + 1] += Nx;
    }
    m_id2X.insert(m_id2X.begin() + id, Nx, -1);
    const int Nd = static_cast<int>(m_id2X.size());
    for (int jd = id + Nx; jd < Nd; ++jd) {
      const int iX = m_id2X[jd];
      if (iX != -1) {
        m_iX2d[iX] += Nx;
      }
    }
    m_id2idx.insert(m_id2idx.begin() + id, Nx, -1);
    m_xs.insert(m_xs.begin() + id, Nx, ::Point2D());
    m_ds.insert(m_ds.begin() + id, Nx, ::Depth::InverseGaussian());
    m_uds.insert(m_uds.begin() + id, Nx, LM_FLAG_TRACK_DEFAULT);
    const int *_ixs = ixs.data() + i1;
    const GlobalMap::Point *Xs = m_IKF.m_Xs.data() + i1;
    ::Point2D *xs = m_xs.data() + id;
    ::Depth::InverseGaussian *ds = m_ds.data() + id;
    for (int i = 0, jd = id; i < Nx; ++i, ++jd) {
      const int ix = _ixs[i], iX = NX1 + ix, idx = KF.Xs[ix].X.idx;
      m_id2X[jd] = iX;
      m_iX2d[iX] = jd;
      m_id2idx[jd] = idx;
      const int _iX = m_idx2iX[idx];
      if (_iX != -1) {
        const int _id = m_iX2d[_iX];
        if (_id != -1) {

#if 1
          m_id2X[_id] = -1;
          m_iX2d[_iX] = -1;
#endif
          m_id2idx[_id] = -1;
        }
      }
      m_idx2iX[idx] = iX;
      const GlobalMap::Point &X = Xs[i];
      xs[i] = X.m_x.m_x;
      ds[i] = X.m_d;
    }
  //}
  return m_IKF;
}

void Internal::ConvertFeatureMeasurements(const std::vector<MapPointMeasurement> &zs,
                                          FRM::Frame *F) {
  F->ClearMeasurements();

  int i, j;
  const float eps = FTR_VARIANCE_EPSILON * m_K.m_K.fxyI();
  const int N = static_cast<int>(m_idx2iX.size()), Nz1 = static_cast<int>(zs.size());
  m_zsSortTmp.resize(Nz1);
  for (i = j = 0; i < Nz1; ++i) {
    const MapPointMeasurement &z1 = zs[i];
    FeatureMeasurement &z2 = m_zsSortTmp[j];
    if (z1.idx >= N) {
      continue;
    }
    const int iX = m_idx2iX[z1.idx];
    if (iX == -1) {
      continue;
    }
    if ((z2.m_id = m_iX2d[iX]) == -1) {
      continue;
    }
    {
      Rectify(m_K.m_K, z1.x, &z2.m_z, &z2.m_W, &m_UM, eps);
    }
    ++j;
  }
  const int Nz2 = j;
  m_zsSortTmp.resize(Nz2);
  std::sort(m_zsSortTmp.begin(), m_zsSortTmp.end());

  FTR::Measurement z;
  for (int i1 = 0, i2 = 0, iKF = 0; i1 < Nz2; i1 = i2) {
    const int id1 = m_zsSortTmp[i1].m_id;
    iKF = static_cast<int>(std::upper_bound(m_iKF2d.begin() + iKF, m_iKF2d.end(), id1) -
                                            m_iKF2d.begin()) - 1;
    const int id0 = m_iKF2d[iKF], id2 = m_iKF2d[iKF + 1];

    i2 = static_cast<int>(std::lower_bound(m_zsSortTmp.begin() + i1, m_zsSortTmp.end(), id2) -
                                           m_zsSortTmp.begin());
    m_zsTmp.resize(0);
    for (i = i1; i < i2; ++i) {

      const FeatureMeasurement &_z = m_zsSortTmp[i];
      z.m_ix = _z.m_id - id0;
      {
        z.m_z = _z.m_z;
        z.m_W = _z.m_W;
        {
          z.m_ix = _z.m_id - id0;
          m_zsTmp.push_back(z);
        }
      }
    }
    F->PushFrameMeasurement(iKF, m_zsTmp);
  }
}

void Internal::AssertConsistency() {
  const int N = static_cast<int>(m_Ts.size()), nKFs = static_cast<int>(m_CsKF.size());
  if (N > 0) {
    for (int i = 1; i < N; ++i) {
      UT_ASSERT(m_Ts[i - 1] < m_Ts[i]);
    }
    if (nKFs > 0) {
      UT_ASSERT(m_CsKF.back().m_iFrm <= m_Ts.front().m_iFrm);
    }
  }
  UT_ASSERT(static_cast<int>(m_iKF2d.size()) == nKFs + 1);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    UT_ASSERT(m_iKF2d[iKF] <= m_iKF2d[iKF + 1]);
  }
  const int Nd = m_iKF2d.back();
  int SNd = 0;
  const int NX = static_cast<int>(m_iX2d.size());
  for (int iX = 0; iX < NX; ++iX) {
    const int id = m_iX2d[iX];
    if (id == -1) {
      continue;
    }
    UT_ASSERT(m_id2X[id] == iX);
    ++SNd;
  }
  //UT_ASSERT(SNd == Nd);
  UT_ASSERT(static_cast<int>(m_id2X.size()) == Nd && static_cast<int>(m_id2idx.size()) == Nd);
  for (int id = 0; id < Nd; ++id) {
    const int iX = m_id2X[id], idx = m_id2idx[id];
    if (iX == -1) {
      UT_ASSERT(idx == -1);
    } else {
      UT_ASSERT(m_iX2d[iX] == id);
      //if (idx != -1) {
      //  UT_ASSERT(m_idx2iX[idx] == iX);
      //}
      UT_ASSERT(idx != -1);
      UT_ASSERT(m_idx2iX[idx] == iX);
    }
  }
  const int idxMax = static_cast<int>(m_idx2iX.size() - 1);
  for (int idx = 0; idx <= idxMax; ++idx) {
    const int iX = m_idx2iX[idx];
    if (iX == -1) {
      continue;
    }
    const int id = m_iX2d[iX];
    if (id != -1) {
      UT_ASSERT(m_id2idx[id] == idx);
    }
  }
  UT_ASSERT(static_cast<int>(m_xs.size()) == Nd);
  UT_ASSERT(static_cast<int>(m_ds.size()) == Nd);
  UT_ASSERT(static_cast<int>(m_uds.size()) == Nd);
#ifdef CFG_SERIAL
  if (m_LBA.m_serial && m_GBA.m_serial) {
    UT_ASSERT(m_LBA.m_CsKF.Size() == nKFs && m_GBA.m_Cs.Size() == nKFs);
    for (int iKF = 0; iKF < nKFs; ++iKF) {
      //m_LBA.m_CsKF[iKF].AssertEqual(m_GBA.m_Cs[iKF], 0, "", -1.0f, -1.0f);
      UT_ASSERT(m_LBA.m_CsKF[iKF] == m_GBA.m_Cs[iKF]);
    }

  }
#endif
}

void LoadParameters(const std::string fileName) {
  if (!UT::FileExists(fileName)) {
    return;
  }
  const Configurator cfgor(fileName.c_str());
  //cfgor.Print();
  ::LoadParameters(cfgor);
}

bool SaveCalibration(const std::string fileName, const Calibration &K) {
  FILE *fp = fopen(fileName.c_str(), "wb");
  if (!fp) {
    return false;
  }
  UT::SaveB(K, fp);
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LoadCalibration(const std::string fileName, Calibration *K) {
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (!fp) {
    return false;
  }
  UT::LoadB(*K, fp);
  fclose(fp);
  UT::PrintLoaded(fileName);
  return true;
}

bool LoadCalibration(const std::string fileName, float T[3][4], Intrinsic *K,
                     float *ba, float *bw/*, float *sa*/) {
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    return false;
  }
  if (ba) {
    ba[0] = ba[1] = ba[2] = 0.0f;
  }
  if (bw) {
    bw[0] = bw[1] = bw[2] = 0.0f;
  }
  {
    int i = 0, j = 0;
    float k[3];
    char line[UT_STRING_WIDTH_MAX];
    memset(K->ds, 0, sizeof(K->ds));
    while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
      const int len = int(std::remove(line, std::remove(line, line + strlen(line), 10), ' ') - line);
      if (len == 0 || line[len - 1] == ':')
        continue;
      line[len] = 0;
      switch (i) {
      case 0: {
        if (j == 3 ||
            sscanf(line, "- [%f, %f, %f, %f]", &T[j][0], &T[j][1], &T[j][2], &T[j][3]) == 4) {
          ++j;
        }
        break; }
      case 1: {
        if (sscanf(line, "- [%f, %f, %f]", &k[0], &k[1], &k[2]) == 3) {
          if (j == 0) {
            K->fx = k[0];
            K->cx = k[2];
          } else if (j == 1) {
            K->fy = k[1];
            K->cy = k[2];
          }
          ++j;
        }
        break; }
      case 2: {
        if (sscanf(line, "- [%f]", &K->ds[j]) == 1)
          ++j;
        break; }
      case 3: {
        if (j == 0 && ba) {
          sscanf(line, "accel_bias: [%f, %f, %f]", &ba[0], &ba[1], &ba[2]);
        }
        if (j == 1 && bw) {
          sscanf(line, "gyro_bias: [%f, %f, %f]", &bw[0], &bw[1], &bw[2]);
        }
        //if (j == 2 && sa) {
        //  sscanf(line, "accel_scale: [%f, %f, %f]", &sa[0], &sa[1], &sa[2]);
        //}
        ++j;
        break; }
      }
      if ((i == 0 && j == 4) || (i == 1 && j == 3) || (i == 2 && j == 8) || (i == 3 && j == 2)) {
        j = 0;
        ++i;
      }
    }
    if (i < 3) {
      return false;
    }
    Rotation3D R;
    R.Set(T[0], T[1], T[2]);
    R.MakeOrthogonal();
    R.Get(T[0], T[1], T[2]);
  }
  UT::PrintLoaded(fileName);
  return true;
}

void PrintCalibration(const Calibration &K) {
  Rigid3D T;
  T.Set(&K.Tu[0][0]);
  T.Print("Tcu = ", false);
  UT::Print("K   = %f %f %f %f\n", K.K.fx, K.K.fy, K.K.cx, K.K.cy);
  UT::Print("d   = %f %f %f %f %f %f %f %f\n", K.K.ds[0], K.K.ds[1], K.K.ds[2], K.K.ds[3],
                                               K.K.ds[4], K.K.ds[5], K.K.ds[6], K.K.ds[7]);
}

bool SaveGroundTruth(const std::string fileName, const std::vector<CameraIMUState> &XsGT) {
  return UT::VectorSaveB(XsGT, fileName.c_str());
}

bool LoadGroundTruth(const std::string fileName, std::vector<CameraIMUState> *XsGT) {
  return UT::VectorLoadB(*XsGT, fileName.c_str());
}

bool SaveGroundTruth(const std::string fileName, const std::vector<Depth> &dsGT) {
  return UT::VectorSaveB(dsGT, fileName.c_str());
}

bool LoadGroundTruth(const std::string fileName, std::vector<Depth> *dsGT) {
  if (!UT::VectorLoadB(*dsGT, fileName.c_str())) {
    return false;
  }
  const int N = static_cast<int>(dsGT->size());
  for (int i = 0; i < N; ++i) {
    Depth &d = dsGT->at(i);
    if (d.d != 0.0f) {
      continue;
    }
    d.d = DEPTH_INITIAL;
    d.s2 = DEPTH_VARIANCE_INITIAL;
  }
  return true;
}

static inline bool LoadIndexes(const std::string fileName, std::vector<int> *idxs) {
  idxs->resize(0);
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    return false;
  }
  int idx;
  while (fscanf(fp, "%d", &idx) == 1) {
    idxs->push_back(idx);
  }
  fclose(fp);
  return true;
}

static inline bool SaveIndexes(const std::string fileName, const std::vector<int> &idxs) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  const int N = static_cast<int>(idxs.size());
  for (int i = 0; i < N; ++i) {
    fprintf(fp, "%d\n", idxs[i]);
  }
  fclose(fp);
  return true; 
}

bool SaveKeyFrames(const std::string fileName, const std::vector<int> &iFrms) {
  return SaveIndexes(fileName, iFrms);
}

bool LoadKeyFrames(const std::string fileName, std::vector<int> *iFrms) {
  return LoadIndexes(fileName, iFrms);
}

bool LoadKeyFrames(const std::string fileName, const std::vector<float> &ts,
                   std::vector<ubyte> *kfs, const float dtMax) {
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    kfs->resize(0);
    return false;
  }
  char line[UT_STRING_WIDTH_MAX];
  int kf, iFrm;
  float t, dt1, dt2;
  std::vector<ubyte> &_kfs = *kfs;
  _kfs.assign(ts.size(), 0);
  std::vector<float>::const_iterator it = ts.begin();
  while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
    if (line[0] == '#') {
      continue;
    }
    sscanf(line, "%f %d", &t, &kf);
    it = std::lower_bound(it, ts.end(), t);
    dt1 = it == ts.begin() ? FLT_MAX : (t - *(it - 1));
    dt2 = it == ts.end() ? FLT_MAX : (*it - t);
    if (dt1 > dtMax && dt2 > dtMax) {
      continue;
    } else if (dt1 < dt2) {
      iFrm = static_cast<int>(it - ts.begin()) - 1;
    } else {
      iFrm = static_cast<int>(it - ts.begin());
    }
    _kfs[iFrm] = static_cast<ubyte>(kf);
  }
  if (!_kfs.empty()) {
    _kfs.front() = 1;
  }
  fclose(fp);
  UT::PrintLoaded(fileName);
  return true;
}

bool SaveMapPoints(const std::string fileName, const std::vector<int> &idxs) {
  return SaveIndexes(fileName, idxs);
}

bool LoadMapPoints(const std::string fileName, std::vector<int> *idxs) {
  return LoadIndexes(fileName, idxs);
}

bool SaveCurrentFrame(const std::string fileName, const CurrentFrame &CF, const KeyFrame &KF) {
  FILE *fp = fopen(fileName.c_str(), "wb");
  if (!fp) {
    return false;
  }
  UT::SaveB(CF.iFrm, fp);
  UT::SaveB(CF.C, fp);
  UT::VectorSaveB(CF.zs, fp);
  UT::VectorSaveB(CF.us, fp);
  UT::SaveB(CF.t, fp);
  UT::SaveB(CF.d, fp);
  UT::StringSaveB(CF.fileName, fp);
  UT::SaveB(KF.iFrm, fp);
  UT::SaveB(KF.C, fp);
  UT::VectorSaveB(KF.zs, fp);
  const int NX = static_cast<int>(KF.Xs.size());
  UT::SaveB(NX, fp);
  for (int iX = 0; iX < NX; ++iX) {
    const MapPoint &X = KF.Xs[iX];
    UT::VectorSaveB(X.zs, fp);
    UT::SaveB(X.X, fp);
  }
  UT::SaveB(KF.d, fp);
  fclose(fp);
  return true;
}

#ifdef IBA_DEBUG_REMOVE_RIGHT_MEASUREMENTS
static inline void RemoveRightMeasurements(std::vector<MapPointMeasurement> *zs) {
  int i, j;
  const int N = static_cast<int>(zs->size());
  for (i = j = 0; i < N; ++i) {
    const MapPointMeasurement &z = zs->at(i);
    if (!z.right) {
      zs->at(j++) = z;
    }
  }
  zs->resize(j);
}
#endif

bool LoadCurrentFrame(const std::string fileName, CurrentFrame *CF, KeyFrame *KF) {
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (!fp) {
    return false;
  }
  UT::LoadB(CF->iFrm, fp);
  UT::LoadB(CF->C, fp);
  UT::VectorLoadB(CF->zs, fp);
  UT::VectorLoadB(CF->us, fp);
  UT::LoadB(CF->t, fp);
  UT::LoadB(CF->d, fp);
  UT::StringLoadB(CF->fileName, fp);
  UT::LoadB(KF->iFrm, fp);
  UT::LoadB(KF->C, fp);
  UT::VectorLoadB(KF->zs, fp);
  const int NX = UT::LoadB<int>(fp);
  KF->Xs.resize(NX);
  for (int iX = 0; iX < NX; ++iX) {
    MapPoint &X = KF->Xs[iX];
    UT::VectorLoadB(X.zs, fp);
    UT::LoadB(X.X, fp);
  }
  UT::LoadB(KF->d, fp);
  fclose(fp);
#ifdef IBA_DEBUG_INVALIDATE_INITIAL_CAMERAS
  CF->C.C.R[0][0] = FLT_MAX;
  KF->C.R[0][0] = FLT_MAX;
#endif
#ifdef IBA_DEBUG_INVALIDATE_INITIAL_POINTS
  for (int iX = 0; iX < NX; ++iX) {
    KF->Xs[iX].X.X[0] = FLT_MAX;
  }
#endif

#ifdef IBA_DEBUG_REMOVE_RIGHT_MEASUREMENTS
  RemoveRightMeasurements(&CF->zs);
  if (KF->iFrm != -1) {
    RemoveRightMeasurements(&KF->zs);
    const int NX = static_cast<int>(KF->Xs.size());
    for (int iX = 0; iX < NX; ++iX) {
      RemoveRightMeasurements(&KF->Xs[iX].zs);
    }
  }
#endif
  return true;
}

bool LoadCurrentFrameTime(const std::string fileName, double *t) {
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (!fp) {
    return false;
  }
  fseek(fp, sizeof(int) + sizeof(CameraIMUState), SEEK_CUR);
  const int Nz = UT::LoadB<int>(fp);
  fseek(fp, sizeof(MapPointMeasurement) * Nz, SEEK_CUR);
  const int Nu = UT::LoadB<int>(fp);
  fseek(fp, sizeof(IMUMeasurement) * Nu, SEEK_CUR);
  *t = static_cast<double>(UT::LoadB<float>(fp));
  fclose(fp);
  return true;
}

bool LoadRelativeConstraints(const std::string fileName, const std::vector<float> &ts,
                             std::vector<RelativeConstraint> *Zs, const float dtMax) {
  Zs->resize(0);
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    return false;
  }
  char line[UT_STRING_WIDTH_MAX];
  float t1, t2;
  RelativeConstraint Z;
  //const float sr = 1.0f * UT_FACTOR_DEG_TO_RAD;
  //const float sr2 = sr * sr;
  std::string format = "%f";
  for (int i = 0; i < 40; ++i) {
    format = format + " %f";
  }
  memset(&Z.S.S[0][0], 0, sizeof(Z.S));
  while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
    if (line[0] == '#') {
      continue;
    }
    sscanf(line, format.c_str(),
           &t1, &t2,
           &Z.T.p[0], &Z.T.p[1], &Z.T.p[2],
           &Z.T.R[0][0], &Z.T.R[0][1], &Z.T.R[0][2],
           &Z.T.R[1][0], &Z.T.R[1][1], &Z.T.R[1][2],
           &Z.T.R[2][0], &Z.T.R[2][1], &Z.T.R[2][2],
           &Z.S.S[0][0], &Z.S.S[0][1], &Z.S.S[0][2],
           &Z.S.S[1][0], &Z.S.S[1][1], &Z.S.S[1][2],
           &Z.S.S[2][0], &Z.S.S[2][1], &Z.S.S[2][2],
           &Z.S.S[0][3], &Z.S.S[0][4], &Z.S.S[0][5],
           &Z.S.S[1][3], &Z.S.S[1][4], &Z.S.S[1][5],
           &Z.S.S[2][3], &Z.S.S[2][4], &Z.S.S[2][5],
           &Z.S.S[3][3], &Z.S.S[3][4], &Z.S.S[3][5],
           &Z.S.S[4][3], &Z.S.S[4][4], &Z.S.S[4][5],
           &Z.S.S[5][3], &Z.S.S[5][4], &Z.S.S[5][5]);
    //Z.S.S[3][3] = Z.S.S[4][4] = Z.S.S[5][5] = sr2;
    for (int i = 0; i < 3; ++i) {
      for (int j = 3; j < 6; ++j) {
        Z.S.S[j][i] = Z.S.S[i][j];
      }
    }
    const std::vector<float>::const_iterator it1 = std::lower_bound(ts.begin(), ts.end(), t1);
    const std::vector<float>::const_iterator it2 = std::lower_bound(ts.begin(), ts.end(), t2);
    if (it1 != ts.end() && fabs(t1 - *it1) <= dtMax && it2 != ts.end() && fabs(t2 - *it2) <= dtMax) {
      Z.iFrm1 = static_cast<int>(it1 - ts.begin());
      Z.iFrm2 = static_cast<int>(it2 - ts.begin());
      Zs->push_back(Z);
    }
  }
  fclose(fp);
  UT::PrintLoaded(fileName);
  return true;
}

void PrintCameraPose(const int iFrm, const CameraPose &C, const bool n) {
  const ::Point3D p = C.p;
  const Quaternion q = Rotation3D(C.R).GetQuaternion();
  if (n) {
    UT::Print("[%d]   %e %e %e   %e %e %e %e\n", iFrm, p.x(), p.y(), p.z(),
              q.x(), q.y(), q.z(), q.w());
  } else {
    UT::Print("\r[%d]   %f %f %f   %f %f %f %f", iFrm, p.x(), p.y(), p.z(),
              q.x(), q.y(), q.z(), q.w());
  }
}

void SaveCameraPose(const int iFrm, const CameraPose &C, FILE *fp) {
  const ::Point3D p = C.p;
  const Quaternion q = Rotation3D(C.R).GetQuaternion();
  fprintf(fp, "[%d]   %e %e %e   %e %e %e %e\n", iFrm, p.x(), p.y(), p.z(),
          q.x(), q.y(), q.z(), q.w());
}

bool SaveCameraPoses(const std::string fileName, const std::vector<int> &iFrms,
                     const std::vector<CameraPose> &Cs) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  const int N = static_cast<int>(iFrms.size());
  for (int i = 0; i < N; ++i) {
    fprintf(fp, "%d\n", iFrms[i]);
    const CameraPose C = Cs[i];
    fprintf(fp, "%e %e %e\n", C.R[0][0], C.R[0][1], C.R[0][2]);
    fprintf(fp, "%e %e %e\n", C.R[1][0], C.R[1][1], C.R[1][2]);
    fprintf(fp, "%e %e %e\n", C.R[2][0], C.R[2][1], C.R[2][2]);
    fprintf(fp, "%e %e %e\n", C.p[0], C.p[1], C.p[2]);
  }
  fclose(fp);
  return true;
}

bool LoadCameraPoses(const std::string fileName, std::vector<int> *iFrms,
                     std::vector<CameraPose> *Cs) {
  iFrms->resize(0);
  Cs->resize(0);
  FILE *fp = fopen(fileName.c_str(), "r");
  if (!fp) {
    return false;
  }
  int iFrm;
  CameraPose C;
  while (fscanf(fp, "%d", &iFrm) == 1 &&
         fscanf(fp, "%e %e %e", &C.R[0][0], &C.R[0][1], &C.R[0][2]) == 3 &&
         fscanf(fp, "%e %e %e", &C.R[1][0], &C.R[1][1], &C.R[1][2]) == 3 &&
         fscanf(fp, "%e %e %e", &C.R[2][0], &C.R[2][1], &C.R[2][2]) == 3 &&
         fscanf(fp, "%e %e %e", &C.p[0], &C.p[1], &C.p[2]) == 3) {
    iFrms->push_back(iFrm);
    Cs->push_back(C);
  }
  fclose(fp);
  return true;
}

void PrintPoints(const std::vector<Point3D> &Xs) {
  const int N = static_cast<int>(Xs.size());
  for (int i = 0; i < N; ++i) {
    const Point3D &X = Xs[i];
    if (X.idx != -1) {
      UT::Print("%d %f %f %f\n", X.idx, X.X[0], X.X[1], X.X[2]);
    }
  }
}

static inline bool AssertError(const int iFrm, const float e, const float eMax, const std::string str) {
  if (e <= eMax) {
    return true;
  }
  UT::Print("[%d] %s = %f > %f\n", iFrm, str.c_str(), e, eMax);
  return false;
}

bool AssertError(const int iFrm, const Error &e) {
  const bool ex = AssertError(iFrm, e.ex, MAX_ERROR_FEATURE, "ex");
  const bool eur = AssertError(iFrm, e.eur, MAX_ERROR_IMU_ROTATION, "eur");
  const bool eup = AssertError(iFrm, e.eup, MAX_ERROR_IMU_POSITION, "eup");
  const bool euv = AssertError(iFrm, e.euv, MAX_ERROR_IMU_VELOCITY, "euv");
  const bool euba = AssertError(iFrm, e.euba, MAX_ERROR_IMU_BIAS_ACCELERATION, "euba");
  const bool eubw = AssertError(iFrm, e.eubw, MAX_ERROR_IMU_BIAS_GYROSCOPE, "eubw");
  const bool edr = AssertError(iFrm, e.edr, MAX_ERROR_DRIFT_ROTATION, "edr");
  const bool edp = AssertError(iFrm, e.edp, MAX_ERROR_DRIFT_POSITION, "edp");
  return ex && eur && eup && euv && euba && eubw && edr && edp;
}

void PrintError(const std::string str, const Error &e) {
  const std::string _str(str.size(), ' ');
  UT::PrintSeparator();
  UT::Print("%sex = %f\n", str.c_str(), e.ex);
  UT::Print("%seu = %f %f %f %f %f\n", _str.c_str(), e.eur, e.eup, e.euv, e.euba, e.eubw);
  UT::Print("%sed = %f %f\n", _str.c_str(), e.edr, e.edp);
}

}  // namespace IBA
