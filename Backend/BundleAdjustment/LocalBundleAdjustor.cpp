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
#include "LocalBundleAdjustor.h"
#include "GlobalBundleAdjustor.h"
#include "IBA_internal.h"
#include "Vector12.h"

#ifdef CFG_DEBUG
#ifdef CFG_DEBUG_EIGEN
#define LBA_DEBUG_EIGEN
#endif
#if WIN32
#define LBA_DEBUG_CHECK
//#define LBA_DEBUG_PRINT
//#define LBA_DEBUG_PRINT_STEP
#ifdef CFG_GROUND_TRUTH
//#define LBA_DEBUG_GROUND_TRUTH_MEASUREMENT
#endif
//#define LBA_DEBUG_PRINT_MARGINALIZATION
#define LBA_DEBUG_VIEW
#ifdef LBA_DEBUG_VIEW
#include "ViewerIBA.h"
static ViewerIBA *g_viewer = NULL;
#endif
#endif
#endif

#ifdef LBA_DEBUG_CHECK
static const int g_ic = 0;
//static const int g_ic = INT_MAX;
static const int g_iKF = 1;
//static const int g_iKF = INT_MAX;
static const int g_ix = 8;
static const int g_ist = 41;
static const int g_izKF = 0;
static const int g_izLF = 22;
static const Rigid3D *g_CKF = NULL;
static const LocalBundleAdjustor::KeyFrame *g_KF = NULL;
static const Depth::InverseGaussian *g_d = NULL;
static const FTR::Factor::FixSource::Source::A *g_Ax = NULL;
static const FTR::Factor::FixSource::Source::A *g_AxST = NULL;
static const FTR::Factor::FixSource::Source::M *g_MxST = NULL;
static const FTR::Factor::FixSource::Source::M *g_Mx = NULL;
static const FTR::Factor::Depth *g_AzKF = NULL;
static const Camera *g_CLF = NULL;
static const LocalBundleAdjustor::LocalFrame *g_LF = NULL;
static const Camera::Factor::Unitary::CC *g_SAcuLF = NULL;
static const Camera::Factor::Unitary::CC *g_SMcuLF = NULL;
static const Camera::Factor *g_SAcmLF = NULL;
static const FTR::Factor::FixSource::A1 *g_Az1LF = NULL;
static const FTR::Factor::FixSource::A2 *g_Az2LF = NULL;
static const FTR::Factor::FixSource::M1 *g_Mz1LF = NULL;
static const FTR::Factor::FixSource::M2 *g_Mz2LF = NULL;
static const FTR::Factor::DD *g_SmddST = NULL;
#endif

void LocalBundleAdjustor::Initialize(IBA::Solver *solver, const int serial, const int verbose,
                                     const int debug, const int history) {
  MT::Thread::Initialize(serial, 4, "LBA");
  m_solver = solver;
  m_LM = &solver->m_internal->m_LM;
  m_GM = &solver->m_internal->m_GM;
  m_GBA = &solver->m_internal->m_GBA;
  m_K = solver->m_internal->m_K;
  m_verbose = verbose;
  m_debug = debug;
  m_history = history;
  m_dir = solver->m_internal->m_dir;
#ifdef CFG_GROUND_TRUTH
  m_CsGT = solver->m_internal->m_CsGT.Data();
#ifdef LBA_DEBUG_GROUND_TRUTH_MEASUREMENT
  if (m_CsGT) {
    LA::AlignedVector3f ba, bw;
    ba.MakeZero();
    bw.MakeZero();
    const int N = solver->m_internal->m_CsGT.Size();
    for (int i = 0; i < N; ++i) {
      const Camera &C = m_CsGT[i];
      ba += C.m_ba;
      bw += C.m_bw;
    }
    const float s = 1.0f / N;
    ba *= s;
    bw *= s;
    Camera *Cs = (Camera *) m_CsGT;
    for (int i = 0; i < N; ++i) {
      Camera &C = Cs[i];
      C.m_ba = ba;
      C.m_bw = bw;
    }
  }
#endif
  m_dsGT = solver->m_internal->m_DsGT.empty() ? NULL : &solver->m_internal->m_dsGT;
#ifdef CFG_HISTORY
  if (m_history >= 3 && (!m_CsGT || !m_dsGT)) {
    m_history = 2;
  }
#endif
#endif
}

void LocalBundleAdjustor::Reset() {
  MT::Thread::Reset();
  MT_WRITE_LOCK_BEGIN(m_MT, MT_TASK_NONE, MT_TASK_LBA_Reset);
  m_ITs1.resize(0);
  m_ITs2.resize(0);
  m_ILFs1.resize(0);
  m_ILFs2.resize(0);
  m_IKFs1.resize(0);
  m_IKFs2.resize(0);
  m_IDKFs1.resize(0);
  m_IDKFs2.resize(0);
  m_IDMPs1.resize(0);
  m_IDMPs2.resize(0);
  m_IUCs1.resize(0);
  m_IUCs2.resize(0);
  MT_WRITE_LOCK_END(m_MT, MT_TASK_NONE, MT_TASK_LBA_Reset);

  m_GM->LBA_Reset();

  for (int i = 0; i < TM_TYPES; ++i) {
    m_ts[i].Reset(TIME_AVERAGING_COUNT);
  }
#ifdef CFG_HISTORY
  m_hists.resize(0);
#endif

  m_delta2 = BA_DL_RADIUS_INITIAL;

  m_ic2LF.reserve(LBA_MAX_LOCAL_FRAMES);    m_ic2LF.resize(0);
  m_LFs.reserve(LBA_MAX_LOCAL_FRAMES);      m_LFs.resize(0);
  m_CsLF.Reserve(LBA_MAX_LOCAL_FRAMES);     m_CsLF.Resize(0);
#ifdef CFG_GROUND_TRUTH
  m_CsLFGT.Reserve(LBA_MAX_LOCAL_FRAMES);   m_CsLFGT.Resize(0);
#endif
  m_ucsLF.reserve(LBA_MAX_LOCAL_FRAMES);    m_ucsLF.resize(0);
  m_ucmsLF.reserve(LBA_MAX_LOCAL_FRAMES);   m_ucmsLF.resize(0);
#ifdef CFG_INCREMENTAL_PCG
  m_xcsLF.Reserve(LBA_MAX_LOCAL_FRAMES);    m_xcsLF.Resize(0);
  m_xmsLF.Reserve(LBA_MAX_LOCAL_FRAMES);    m_xmsLF.Resize(0);
#endif
  m_DsLF.Reserve(LBA_MAX_LOCAL_FRAMES);     m_DsLF.Resize(0);
#ifdef CFG_GROUND_TRUTH
  m_DsLFGT.Reserve(LBA_MAX_LOCAL_FRAMES);   m_DsLFGT.Resize(0);
#endif
  m_AdsLF.Reserve(LBA_MAX_LOCAL_FRAMES);    m_AdsLF.Resize(0);
  m_AfpsLF.Reserve(LBA_MAX_LOCAL_FRAMES);   m_AfpsLF.Resize(0);
  m_AfmsLF.Reserve(LBA_MAX_LOCAL_FRAMES);   m_AfmsLF.Resize(0);
  m_SAcusLF.Reserve(LBA_MAX_LOCAL_FRAMES);  m_SAcusLF.Resize(0);
  m_SMcusLF.Reserve(LBA_MAX_LOCAL_FRAMES);  m_SMcusLF.Resize(0);
  m_SAcmsLF.Reserve(LBA_MAX_LOCAL_FRAMES);  m_SAcmsLF.Resize(0);

  m_KFs.resize(0);
  m_iFrmsKF.resize(0);
  m_CsKF.Resize(0);
#ifdef CFG_GROUND_TRUTH
  m_CsKFGT.Resize(0);
#endif
  m_ucsKF.resize(0);
#ifdef CFG_HANDLE_SCALE_JUMP
  m_dsKF.resize(0);
#endif
  m_usKF.Resize(0);
  m_usKFLast.Resize(0);

  m_iKF2d.assign(1, 0);
  m_ds.resize(0);
  m_uds.resize(0);

#ifdef CFG_CHECK_REPROJECTION
  m_esLF.resize(0);
  m_esKF.resize(0);
#endif

  m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL, BA_VARIANCE_PRIOR_VELOCITY_FIRST,
                    BA_VARIANCE_PRIOR_BIAS_ACCELERATION_FIRST,
                    BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_FIRST);
#ifdef LBA_DEBUG_GROUND_TRUTH_MEASUREMENT
  if (m_CsGT) {
    m_ZpLF.DebugSetMeasurement(m_CsGT[0]);
  }
#endif
  m_Zp.Initialize(m_ZpLF);
  m_ApLF.MakeZero();
  m_ZpKF.Invalidate();
  //m_F = 0.0f;
}

void LocalBundleAdjustor::PushLocalFrame(const InputLocalFrame &ILF) {
  MT_WRITE_LOCK_BEGIN(m_MT, ILF.m_T.m_iFrm, MT_TASK_LBA_PushLocalFrame);
  m_ITs1.push_back(IT_LOCAL_FRAME);
  m_ILFs1.push_back(ILF);
  MT_WRITE_LOCK_END(m_MT, ILF.m_T.m_iFrm, MT_TASK_LBA_PushLocalFrame);
}

void LocalBundleAdjustor::PushKeyFrame(const GlobalMap::InputKeyFrame &IKF, const bool serial) {
  MT_WRITE_LOCK_BEGIN(m_MT, IKF.m_T.m_iFrm, MT_TASK_LBA_PushKeyFrame);
  if (m_serial) {
    m_ITs1.push_back(IT_KEY_FRAME_SERIAL);
  } else {
    m_ITs1.push_back(IT_KEY_FRAME);
  }
  m_IKFs1.push_back(IKF);
  MT_WRITE_LOCK_END(m_MT, IKF.m_T.m_iFrm, MT_TASK_LBA_PushKeyFrame);
}

void LocalBundleAdjustor::PushDeleteKeyFrame(const int iFrm, const int iKF){
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_PushDeleteKeyFrame);
  m_ITs1.push_back(IT_DELETE_KEY_FRAME);
  m_IDKFs1.push_back(iKF);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_LBA_PushDeleteKeyFrame);
}

void LocalBundleAdjustor::PushDeleteMapPoints(const int iFrm, const std::vector<int> &ids) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_PushDeleteMapPoints);
  m_ITs1.push_back(IT_DELETE_MAP_POINTS);
  m_IDMPs1.push_back(ids);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_LBA_PushDeleteMapPoints);
}

void LocalBundleAdjustor::PushUpdateCameras(const int iFrm,
                                            const std::vector<GlobalMap::InputCamera> &Cs,
                                            const bool serial) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_PushUpdateCameras);
  if (m_serial) {
    m_ITs1.push_back(IT_UPDATE_CAMERAS_SERIAL);
  } else {
    m_ITs1.push_back(IT_UPDATE_CAMERAS);
  }
  m_IUCs1.push_back(Cs);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_LBA_PushUpdateCameras);
}

void LocalBundleAdjustor::GetCamera(FRM::Tag &T, Camera &C) {
  MT_READ_LOCK_BEGIN(m_MTC, MT_TASK_NONE, MT_TASK_LBA_GetCamera);
  T = m_C.m_T;
  C = m_C.m_C;
  MT_READ_LOCK_END(m_MTC, MT_TASK_NONE, MT_TASK_LBA_GetCamera);
}

void LocalBundleAdjustor::Run() {
  m_delta2 = BA_DL_RADIUS_INITIAL;
  m_ts[TM_TOTAL].Start();
  m_ts[TM_SYNCHRONIZE].Start();
  SynchronizeData();
  m_ts[TM_SYNCHRONIZE].Stop();
  m_ts[TM_TOTAL].Stop();

  const int iFrm = m_LFs[m_ic2LF.back()].m_T.m_iFrm;

  m_ts[TM_TOTAL].Start();
  for (m_iIter = 0; m_iIter < BA_MAX_ITERATIONS; ++m_iIter) 
  {
    m_ts[TM_FACTOR].Start();
    UpdateFactors();
    m_ts[TM_FACTOR].Stop();

    m_ts[TM_SCHUR_COMPLEMENT].Start();
    UpdateSchurComplement();
    m_ts[TM_SCHUR_COMPLEMENT].Stop();

    m_ts[TM_CAMERA].Start();
    const bool scc = SolveSchurComplement();
    m_ts[TM_CAMERA].Stop();

    m_ts[TM_DEPTH].Start();
    SolveBackSubstitution();
    m_ts[TM_DEPTH].Stop();

    m_xsGD.Resize(0);   m_x2GD = 0.0f;
    m_xsDL.Resize(0);   m_x2DL = 0.0f;
    m_rho = FLT_MAX;
    const int nItersDL = std::max(BA_DL_MAX_ITERATIONS, 1);
    for (m_iIterDL = 0; m_iIterDL < nItersDL; ++m_iIterDL) {
      if (m_x2GN > m_delta2 && m_xsGD.Empty() && BA_DL_MAX_ITERATIONS > 0) {
        m_ts[TM_UPDATE].Start();
        SolveGradientDescent();
        m_ts[TM_UPDATE].Stop();
      }

      m_ts[TM_UPDATE].Start();
      SolveDogLeg();
      UpdateStatesPropose();
      m_ts[TM_UPDATE].Stop();

      if (BA_DL_MAX_ITERATIONS == 0) {
        break;
      }

      m_ts[TM_UPDATE].Start();
      ComputeReduction();
      m_ts[TM_UPDATE].Stop();

      m_ts[TM_UPDATE].Start();
      const bool accept = UpdateStatesDecide();
      m_ts[TM_UPDATE].Stop();

      if (accept) {
        break;
      }
    }
    if (m_iIterDL == BA_DL_MAX_ITERATIONS) {
      m_ts[TM_UPDATE].Start();
      UpdateStatesDecide();
      m_ts[TM_UPDATE].Stop();
    }

    if (!m_update || m_converge) {
      break;
    }
    if (LBA_EMBEDDED_POINT_ITERATION) {
      m_ts[TM_UPDATE].Start();
      EmbeddedPointIteration(m_CsLF, m_CsKF, m_ucsKF, m_uds, &m_ds);
      m_ts[TM_UPDATE].Stop();
    }
    MT_READ_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_BufferDataEmpty);
    m_empty = BufferDataEmpty();
    MT_READ_LOCK_END(m_MT, iFrm, MT_TASK_LBA_BufferDataEmpty);
    if (!m_empty) {
      break;
    }
  }

  m_ts[TM_TOTAL].Stop();
  for (int i = 0; i < TM_TYPES; ++i) {
    m_ts[i].Finish();
  }


  UpdateData();

  // Trigger LBA callback function if set
  if (m_callback) {
    m_callback(iFrm, m_LFs[m_ic2LF.back()].m_T.m_t);
  }
}

void LocalBundleAdjustor::SetCallback(const IBA::Solver::IbaCallback& iba_callback) {
  m_callback = iba_callback;
}

//int LocalBundleAdjustor::GetTotalPoints(int *N) {
//  int SN = 0;
//  const int _N = static_cast<int>(m_hists.size());
//  for (int i = 0; i < _N; ++i) {
//    SN += m_hists[i].m_Nd;
//  }
//  if (N) {
//    *N = _N;
//  }
//  return SN;
//}

float LocalBundleAdjustor::GetTotalTime(int *N) {
#ifdef CFG_HISTORY
  const int _N = static_cast<int>(m_hists.size());
  m_work.Resize(_N);
  LA::AlignedVectorXf ts(m_work.Data(), _N, false);
  for (int i = 0; i < _N; ++i) {
    ts[i] = static_cast<float>(m_hists[i].m_ts[TM_TOTAL]);
  }
  if (N) {
    *N = _N;
  }
  return ts.Sum();
#else
  if (N) {
    *N = 0;
  }
  return 0.0f;
#endif
}

bool LocalBundleAdjustor::SaveTimes(const std::string fileName) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const double *ts = m_hists[i].m_ts;
    for (int j = 0; j < TM_TYPES; ++j) {
      fprintf(fp, "%f ", ts[j]);
    }
    fprintf(fp, "\n");
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LocalBundleAdjustor::SaveCameras(const std::string fileName, const bool poseOnly) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  Point3D p;
  Quaternion q;
  Rotation3D R;
  LA::AlignedVector3f ba, bw;
  const Rotation3D RuT = m_K.m_Ru.GetTranspose();
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const History &hist = m_hists[i];
    const Rigid3D &T = hist.m_C.m_T;
    p = hist.m_C.m_p + T.GetAppliedRotationInversely(m_K.m_pu);
    Rotation3D::AB(RuT, T, R);
    R.GetQuaternion(q);
    fprintf(fp, "%f %f %f %f %f %f %f %f", hist.m_t, p.x(), p.y(), p.z(),
                                           q.x(), q.y(), q.z(), q.w());
    if (!poseOnly) {
      RuT.Apply(hist.m_C.m_ba, ba);
      RuT.Apply(hist.m_C.m_bw, bw);
      ba += m_K.m_ba;
      bw += m_K.m_bw;
      const LA::AlignedVector3f &v = hist.m_C.m_v;
      fprintf(fp, " %f %f %f %f %f %f %f %f %f", v.x(), v.y(), v.z(),
                                                 ba.x(), ba.y(), ba.z(),
                                                 bw.x(), bw.y(), bw.z());
    }
    fprintf(fp, "\n");
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LocalBundleAdjustor::SaveCosts(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const History &hist = m_hists[i];
    switch (type) {
    case 0: hist.m_ESa.Save(fp);    break;
    case 1: hist.m_ESb.Save(fp);    break;
    case 2: hist.m_ESp.Save(fp);    break;
#ifdef CFG_GROUND_TRUTH
    case 3: hist.m_ESaGT.Save(fp);  break;
    case 4: hist.m_ESpGT.Save(fp);  break;
#endif
    }
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LocalBundleAdjustor::SaveResiduals(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const History &hist = m_hists[i];
    switch (type) {
    case 0: hist.m_R.Save(fp);    break;
#ifdef CFG_GROUND_TRUTH
    case 1: hist.m_RGT.Save(fp);  break;
#endif
    }
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LocalBundleAdjustor::SavePriors(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const History &hist = m_hists[i];
    switch (type) {
#ifdef CFG_GROUND_TRUTH
    case 0: hist.m_PS.Save(fp);   break;
    case 1: hist.m_PSKF.Save(fp); break;
    case 2: hist.m_PSLF.Save(fp); break;
#endif
    }
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool LocalBundleAdjustor::SaveMarginalizations(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
#ifdef CFG_HISTORY
  const int N = static_cast<int>(m_hists.size());
  for (int i = 0; i < N; ++i) {
    const History &hist = m_hists[i];
    switch (type) {
    case 0: hist.m_MSLP.Save(fp); break;
    case 1: hist.m_MSEM.Save(fp); break;
#ifdef CFG_GROUND_TRUTH
    case 2: hist.m_MSGT.Save(fp); break;
#endif
    }
  }
#endif
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

void LocalBundleAdjustor::ComputeErrorFeature(float *ex) {
  float exi;
  *ex = 0.0f;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    ComputeErrorFeature(&m_LFs[iLF], m_CsLF[iLF].m_T, m_CsKF, m_ds, &exi);
    *ex = std::max(exi, *ex);
  }
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    ComputeErrorFeature(&m_KFs[iKF], m_CsKF[iKF], m_CsKF, m_ds, &exi, iKF);
    *ex = std::max(exi, *ex);
  }
}

void LocalBundleAdjustor::ComputeErrorFeature(const FRM::Frame *F, const Rigid3D &C,
                                              const AlignedVector<Rigid3D> &CsKF,
                                              const std::vector<Depth::InverseGaussian> &ds,
                                              float *ex, const int iKF) {
  Rigid3D Tr[2];
  FTR::Error e;
//#if 0
#if 1
  float Se2 = 0.0f;
#ifdef CFG_STEREO
  float Se2r = 0.0f;
#endif
  int SN = 0;
#else
  const int Nz = static_cast<int>(F->m_zs.size());
#ifdef CFG_STEREO
  m_work.Resize(Nz * 2);
  LA::AlignedVectorXf e2s(m_work.Data(), Nz, false);
  LA::AlignedVectorXf e2rs(e2s.BindNext(), Nz, false);
  e2s.Resize(0);
  e2rs.Resize(0);
#else
  m_work.Resize(Nz);
  LA::AlignedVectorXf e2s(m_work.Data(), Nz, false);
  e2s.Resize(0);
#endif
#endif
  const int NZ = static_cast<int>(F->m_Zs.size());
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = F->m_Zs[iZ];
    *Tr = C / CsKF[Z.m_iKF];
#ifdef CFG_STEREO
    Tr[1] = Tr[0];
    Tr[1].SetTranslation(m_K.m_br + Tr[0].GetTranslation());
#endif
    const Depth::InverseGaussian *_ds = ds.data() + m_iKF2d[Z.m_iKF];
    const KeyFrame &KF = m_KFs[Z.m_iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      const FTR::Measurement &z = F->m_zs[iz];
      const int ix = z.m_ix;
      FTR::GetError(Tr, KF.m_xs[ix], _ds[ix], z, e);
#ifdef CFG_STEREO
      if (z.m_z.Valid()) {
        const float e2 = e.m_e.SquaredLength();
#if 1
//#if 0
        Se2 = e2 + Se2;
        ++SN;
#else
        e2s.Push(e2);
#endif
      }
      if (z.m_zr.Valid()) {
        const float e2r = e.m_er.SquaredLength();
#if 1
//#if 0
        Se2r = e2r + Se2r;
        ++SN;
#else
        e2rs.Push(e2r);
#endif
      }
#else
      const float e2 = e.m_e.SquaredLength();
#if 1
//#if 0
      Se2 = e2 * m_K.m_K.fxy() + Se2;
      ++SN;
#else
      e2s.Push(e2);
#endif
#endif
    }
  }
#ifdef CFG_STEREO
  if (iKF != -1) {
    const Depth::InverseGaussian *_ds = ds.data() + m_iKF2d[iKF];
    const KeyFrame &KF = *((KeyFrame *) F);
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (KF.m_xs[ix].m_xr.Invalid()) {
        continue;
      }
      FTR::GetError(m_K.m_br, _ds[ix], KF.m_xs[ix], e.m_er);
      const float e2r = e.m_er.SquaredLength();
#if 1
//#if 0
      Se2r = e2r + Se2r;
      ++SN;
#else
      e2rs.Push(e2r);
#endif
    }
  }
#endif
//#if 0
#if 1
#ifdef CFG_STEREO
  *ex = SN == 0 ? 0.0f : sqrtf((Se2 * m_K.m_K.fxy() + Se2r * m_K.m_Kr.fxy()) / SN);
#else
  *ex = SN == 0 ? 0.0f : sqrtf(Se2 * m_K.m_K.fxy() / SN);
#endif
#else
  if (!e2s.Empty()) {
    const int ith = e2s.Size() >> 1;
    std::nth_element(e2s.Data(), e2s.Data() + ith, e2s.End());
    *ex = sqrtf(e2s[ith] * m_K.m_K.fxy());
  } else
    *ex = 0.0f;
#ifdef CFG_STEREO
  if (!e2rs.Empty()) {
    const int ith = e2rs.Size() >> 1;
    std::nth_element(e2rs.Data(), e2rs.Data() + ith, e2rs.End());
    *ex = std::max(sqrtf(e2rs[ith] * m_K.m_K.fxy()), *ex);
  }
#endif
#endif
}

void LocalBundleAdjustor::ComputeErrorIMU(float *er, float *ep, float *ev,
                                          float *eba, float *ebw) {
  *er = *ep = *ev = *eba = *ebw = 0.0f;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int ic1 = 0, ic2 = 1; ic2 < nLFs; ic1 = ic2++) {
    const int iLF1 = m_ic2LF[ic1], iLF2 = m_ic2LF[ic2];
    const IMU::Delta::Error e = m_DsLF[iLF2].GetError(m_CsLF[iLF1], m_CsLF[iLF2], m_K.m_pu,
                                                      BA_ANGLE_EPSILON);
    *er = std::max(e.m_er.SquaredLength(), *er);
    *ep = std::max(e.m_ep.SquaredLength(), *ep);
    *ev = std::max(e.m_ev.SquaredLength(), *ev);
    *eba = std::max(e.m_eba.SquaredLength(), *eba);
    *ebw = std::max(e.m_ebw.SquaredLength(), *ebw);
  }
  *er *= UT_FACTOR_RAD_TO_DEG;
  *ebw *= UT_FACTOR_RAD_TO_DEG;
}

void LocalBundleAdjustor::ComputeErrorDrift(float *er, float *ep) {
  *er = 0.0f;
  *ep = 0.0f;
}

float LocalBundleAdjustor::ComputeRMSE() {
  float Se2 = 0.0f;
  //Changed by H.K: ‘m_hists’ was not declared in this scope
  const int N = 1;//static_cast<int>(m_hists.size());
  return sqrtf(Se2 / N);
}

float LocalBundleAdjustor::GetTotalDistance() {
  float S = 0.0f;
  return S;
}

void LocalBundleAdjustor::SynchronizeData() {
  const int iFrm = m_ic2LF.empty() ? MT_TASK_NONE : m_LFs[m_ic2LF.back()].m_T.m_iFrm;
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_SynchronizeData);
  m_ITs1.swap(m_ITs2);
  m_ILFs1.swap(m_ILFs2);
  m_IKFs1.swap(m_IKFs2);
  m_IDKFs1.swap(m_IDKFs2);
  m_IDMPs1.swap(m_IDMPs2);
  m_IUCs1.swap(m_IUCs2);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_LBA_SynchronizeData);
  m_UcsLF.assign(m_LFs.size(), LM_FLAG_FRAME_DEFAULT);
  m_UcsKF.assign(m_KFs.size(), LM_FLAG_FRAME_DEFAULT);
  m_Uds.assign(m_ds.size(), LM_FLAG_TRACK_DEFAULT);

  for (std::list<InputLocalFrame>::iterator ILF = m_ILFs2.begin(); ILF != m_ILFs2.end(); ++ILF) {
    if (ILF->m_C.m_T.Valid() && ILF->m_C.m_v.Valid()) {
      continue;
    }
    Camera C;
    if (m_LFs.empty()) {
      IMU::InitializeCamera(ILF->m_us, C);
    } else {
      float _t;
      Camera _C;
      const IMU::Measurement *_u;
      if (ILF == m_ILFs2.begin()) {
        const int iLF = m_ic2LF.back();
        const LocalFrame &LF = m_LFs[iLF];
        _t = LF.m_T.m_t;
        _C = m_CsLF[iLF];
        _u = &LF.m_us.Back();
      } else {
        std::list<InputLocalFrame>::const_iterator _ILF = ILF;
        --_ILF;
        _t = _ILF->m_T.m_t;
        _C = _ILF->m_C;
        _u = &_ILF->m_us.Back();
      }
      if (LBA_PROPAGATE_CAMERA) {
        IMU::Delta D;
        IMU::PreIntegrate(ILF->m_us, _t, ILF->m_T.m_t, _C, &D, &m_work,
                          false, _u, NULL, BA_ANGLE_EPSILON);
        IMU::Propagate(m_K.m_pu, D, _C, C, BA_ANGLE_EPSILON);
      } else {
        C = _C;
      }
    }
    if (ILF->m_C.m_T.Invalid()) {
      ILF->m_C.m_T = C.m_T;
      ILF->m_C.m_p = C.m_p;
    }
    if (ILF->m_C.m_v.Invalid()) {
      ILF->m_C.m_v = C.m_v;
      ILF->m_C.m_ba = C.m_ba;
      ILF->m_C.m_bw = C.m_bw;
    }
  }
  const bool newKF = !m_IKFs2.empty(), delKF = !m_IDKFs2.empty(), updCams = !m_IUCs2.empty();
  bool serialGBA = false;
  if (newKF) {
    const float w = 1.0f;
    //const float w = BA_WEIGHT_FEATURE;
    m_CsKFBkp.Set(m_CsKF);
    std::list<int>::const_iterator IDK = m_IDKFs2.begin();
    std::list<GlobalMap::InputKeyFrame>::iterator IKF = m_IKFs2.begin();
    for (std::list<InputType>::const_iterator IT = m_ITs2.begin(); IT != m_ITs2.end(); ++IT) {
      if (*IT == IT_DELETE_KEY_FRAME) {
        m_CsKFBkp.Erase(*IDK);
        ++IDK;
      }
      if (*IT != IT_KEY_FRAME && *IT != IT_KEY_FRAME_SERIAL) {
        continue;
      }
      LA::Vector3f Rx;
      LA::SymmetricMatrix2x2f W;
      const int NX = static_cast<int>(IKF->m_Xs.size());
      if (IKF->m_C.Invalid()) {
        for (std::list<InputLocalFrame>::iterator ILF = m_ILFs2.begin(); ILF != m_ILFs2.end(); ++ILF) {
          if (ILF->m_T == IKF->m_T) {
            IKF->m_C = ILF->m_C;
            break;
          }
        }
      }
      if (IKF->m_C.Invalid()) {
        const int nLFs = static_cast<int>(m_LFs.size());
        for (int iLF = 0; iLF < nLFs; ++iLF) {
          if (m_LFs[iLF].m_T == IKF->m_T) {
            IKF->m_C = m_CsLF[iLF];
          }
        }
      }
      m_CsKFBkp.Push(IKF->m_C.m_T);
      const int nKFs = m_CsKFBkp.Size();
      m_R12s.Resize(nKFs);
      m_t12s.Resize(nKFs);

      for (int iX1 = 0, iX2 = 0; iX1 < NX; iX1 = iX2) {
        const int iKF = IKF->m_Xs[iX1].m_iKF;
        const Rigid3D &C = m_CsKFBkp[iKF];
        m_marksTmp1.assign(nKFs, 0);
        for (iX2 = iX1 + 1; iX2 < NX && IKF->m_Xs[iX2].m_iKF == iKF; ++iX2) {}
        for (int iX = iX1; iX < iX2; ++iX) {
          GlobalMap::Point &X = IKF->m_Xs[iX];
          m_zds.resize(0);

          const int Nz = static_cast<int>(X.m_zs.size());
          for (int i = 0; i < Nz; ++i) {
            const FTR::Measurement &z = X.m_zs[i];
            Rotation3D &R = m_R12s[z.m_iKF];

            LA::AlignedVector3f *t = m_t12s.Data() + z.m_iKF;

            if (!m_marksTmp1[z.m_iKF]) {
              const Rigid3D &_C = m_CsKFBkp[z.m_iKF];
              if (_C.Valid()) {
                const Rigid3D T = _C / C;
                R = T;
                T.GetTranslation(*t);
              } else {
                R.Invalidate();
                t->Invalidate();
              }
              m_marksTmp1[z.m_iKF] = 1;
            }
            if (R.Valid()) {
              R.Apply(X.m_x.m_x, Rx);
              z.m_W.GetScaled(w, W);
              m_zds.push_back(Depth::Measurement(*t, Rx, z.m_z, W));

            }
          }

          float eAvg;
          const int Nzd = static_cast<int>(m_zds.size());
          if (!Depth::Triangulate(w, Nzd, m_zds.data(), &X.m_d, &m_work, X.m_d.Valid(), &eAvg) ||
              m_K.m_K.fx() * eAvg > DEPTH_TRI_MAX_ERROR) {
            if (Nzd == 0) {
              X.m_d.Invalidate();
            } else {
              X.m_d.Initialize();
              for (int n = 1; n <= Nzd; ++n) {
                Depth::Triangulate(w, n, m_zds.data(), &X.m_d, &m_work);
              }
              if (!X.m_d.Valid()) {
                X.m_d.Invalidate();
              }
            }
          }
        }
      }
      ++IKF;
    }
  }
  while (!m_ITs2.empty()) {
    const int nKFs = static_cast<int>(m_KFs.size());
    const InputType IT = m_ITs2.front();
    m_ITs2.pop_front();
    if (IT == IT_LOCAL_FRAME) {
      InputLocalFrame &ILF = m_ILFs2.front();
      if (ILF.m_d.Invalid()) {
        if (m_LFs.empty()) {
          ILF.m_d.Initialize();
        } else {
          const LocalFrame &LF = m_LFs[m_ic2LF.back()];
          ILF.m_d = LF.m_d;
          ILF.m_d.Propagate(ILF.m_T.m_t - LF.m_T.m_t);
        }
        const int Nz = static_cast<int>(ILF.m_zs.size()), NzC = SIMD_FLOAT_CEIL(Nz);
        m_work.Resize(NzC * 3);
        LA::AlignedVectorXf us(m_work.Data(), Nz, false);     us.Resize(0);
        LA::AlignedVectorXf ws(us.Data() + NzC, Nz, false);   ws.Resize(0);
        LA::AlignedVectorXf wus(ws.Data() + NzC, Nz, false);
        for (std::list<GlobalMap::InputKeyFrame>::iterator IKF = m_IKFs2.begin();
             IKF != m_IKFs2.end(); ++IKF) {
          if (IKF->m_T != ILF.m_T) {
            continue;
          }
          const int NX = static_cast<int>(IKF->m_Xs.size()), N = NX + Nz, NC = SIMD_FLOAT_CEIL(N);
          m_work.Resize(NC * 3);
          us.Bind(m_work.Data(), N);     us.Resize(0);
          ws.Bind(us.Data() + NC, N);    ws.Resize(0);
          wus.Bind(ws.Data() + NC, N);
          for (int iX = 0; iX < NX; ++iX) {
            const GlobalMap::Point &X = IKF->m_Xs[iX];
            if (X.m_iKF != nKFs || !X.m_d.Valid()) {
              continue;
            }
            us.Push(X.m_d.u());
            ws.Push(X.m_d.s2());
          }
          break;
        }
        Depth::InverseGaussian dz;
        Rigid3D::Row Crz;
        const Rigid3D C = ILF.m_C.m_T;
        const Rigid3D::Row Cz = C.GetRowZ();
        const int NZ = static_cast<int>(ILF.m_Zs.size());
        for (int iZ = 0; iZ < NZ; ++iZ) {
          const FRM::Measurement &Z = ILF.m_Zs[iZ];
          if (Z.m_iKF == nKFs) {
            continue;
          }
          const Depth::InverseGaussian *ds = m_ds.data() + m_iKF2d[Z.m_iKF];
          Rigid3D::ABI(Cz, m_CsKF[Z.m_iKF], Crz);
          const GlobalMap::KeyFrame &KF = m_KFs[Z.m_iKF];
          for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
            const int ix = ILF.m_zs[iz].m_ix;
            const Depth::InverseGaussian &d = ds[ix];
            if (!d.Valid() || !d.ProjectD(Crz, KF.m_xs[ix].m_x, dz)) {
              continue;
            }
            us.Push(dz.u());
            ws.Push(dz.s2());
          }
        }
        ws += DEPTH_VARIANCE_EPSILON;
        ws.MakeInverse();
        const float Sw = ws.Sum();
        if (Sw < FLT_EPSILON) {
          ws.Set(1.0f / Nz);
        } else {
          ws *= 1.0f / Sw;
        }
        LA::AlignedVectorXf::AB(ws, us, wus);
        const float u = wus.Sum();
        us -= u;
        us.MakeSquared();
        LA::AlignedVectorXf::AB(ws, us, wus);
        const float s2 = wus.Sum();
        if (s2 > 0.0f) {
          ILF.m_d.Update(Depth::InverseGaussian(u, s2));
        }
      }
      //if (!ILF.m_Zs.empty() && ILF.m_iKFsMatch.empty()) {
        SearchMatchingKeyFrames(ILF);
      //}
      const int Nk = static_cast<int>(ILF.m_iKFsMatch.size());
      if (ILF.m_iKFNearest == -1 && Nk != 0) {
        if (LBA_MARGINALIZATION_REFERENCE_NEAREST) {
          const int iKFMax = ILF.m_iKFsMatch.back();
          ubyte first = 1;
          int iKFNearest = m_Zp.m_iKFr == -1 ? nKFs - 1 : m_Zp.m_iKFr;
          float imgMotionNearest = FLT_MAX;
          const Rigid3D C = ILF.m_C.m_T;
          const float z = 1.0f / ILF.m_d.u();
          m_marksTmp1.assign(iKFMax + 1, 0);
          for (int ik = 0; ik < Nk; ++ik) {
            m_marksTmp1[ILF.m_iKFsMatch[ik]] = 1;
          }
          const int _nKFs = std::min(nKFs, iKFMax + 1);
          for (int iKF = 0; iKF < _nKFs; ++iKF) {
            if (!m_marksTmp1[iKF]) {
              continue;
            }
            const Rigid3D _C = m_CsKF[iKF];
            const float imgMotion = ComputeImageMotion(z, C, _C, &first);
            if (imgMotion > imgMotionNearest) {
              continue;
            }
            imgMotionNearest = imgMotion;
            iKFNearest = iKF;
          }
          ILF.m_iKFNearest = iKFNearest;
        } else {
          ILF.m_iKFNearest = nKFs - 1;
        }
      } else if (ILF.m_iKFNearest == -1 && Nk == 0) {
        ILF.m_iKFNearest = m_Zp.m_iKFr == -1 ? nKFs - 1 : m_Zp.m_iKFr;
      }
      _PushLocalFrame(ILF);
      m_ILFs2.pop_front();
    } else if (IT == IT_KEY_FRAME || IT == IT_KEY_FRAME_SERIAL) {
      if (IT == IT_KEY_FRAME_SERIAL) {
        serialGBA = true;
      }
      GlobalMap::InputKeyFrame &IKF = m_IKFs2.front();
      const bool v1 = IKF.m_C.m_T.Valid(), v2 = IKF.m_C.m_v.Valid();
      if (!v1 || !v2) {
        const int nLFs = static_cast<int>(m_LFs.size());
        for (int ic = nLFs - 1; ic >= 0; --ic) {
          const int iLF = m_ic2LF[ic];
          if (m_LFs[iLF].m_T != IKF.m_T) {
            continue;
          }
          const Camera &C = m_CsLF[iLF];
          if (!v1) {
            IKF.m_C.m_T = C.m_T;
            IKF.m_C.m_p = C.m_p;
          }
          if (!v2) {
            IKF.m_C.m_v = C.m_v;
            IKF.m_C.m_ba = C.m_ba;
            IKF.m_C.m_bw = C.m_bw;
          }
          break;
        }
      }
      if (IKF.m_d.Invalid()) {
        const int nLFs = static_cast<int>(m_LFs.size());
        for (int ic = nLFs - 1; ic >= 0; --ic) {
          const int iLF = m_ic2LF[ic];
          const LocalFrame &LF = m_LFs[iLF];
          if (LF.m_T != IKF.m_T) {
            continue;
          }
          IKF.m_d = LF.m_d;
          break;
        }
      }
      SearchMatchingKeyFrames(IKF);
      
      const int NX = static_cast<int>(IKF.m_Xs.size());
      for (int iX = 0; iX < NX; ++iX) {
        GlobalMap::Point &X = IKF.m_Xs[iX];
        if (X.m_d.Invalid()) {
          X.m_d.Initialize(X.m_iKF == nKFs ? IKF.m_d.u() : m_KFs[X.m_iKF].m_d.u());
        }
      }
      _PushKeyFrame(IKF);
      m_IKFs2.pop_front();
    } else if (IT == IT_DELETE_KEY_FRAME) {
      DeleteKeyFrame(m_IDKFs2.front());
      m_IDKFs2.pop_front();
    } else if (IT == IT_DELETE_MAP_POINTS) {
      DeleteMapPoints(m_IDMPs2.front());
      m_IDMPs2.pop_front();
    } else if (IT == IT_UPDATE_CAMERAS || IT == IT_UPDATE_CAMERAS_SERIAL) {
      if (IT == IT_UPDATE_CAMERAS_SERIAL) {
        serialGBA = true;
      }
      UpdateCameras(m_IUCs2.front());
      m_IUCs2.pop_front();
    }
  }
  if (newKF || delKF || updCams) {
    m_ts[TM_TOTAL].Stop();
    m_ts[TM_SYNCHRONIZE].Stop();
    //Wake up GBA
    m_GBA->WakeUp(serialGBA);
    m_ts[TM_SYNCHRONIZE].Start();
    m_ts[TM_TOTAL].Start();
  }
  if (!m_iFrmsKF.empty() &&
      m_GM->LBA_Synchronize(m_iFrmsKF.back(), m_CsKF, m_CsKFBkp, m_marksTmp1)) {
    UpdateCameras(m_marksTmp1, m_CsKFBkp, m_CsKF);
  }
}

void LocalBundleAdjustor::UpdateData() {
  const int iFrm1 = m_LFs[m_ic2LF.front()].m_T.m_iFrm;
  const int iFrm2 = m_LFs[m_ic2LF.back()].m_T.m_iFrm;
#ifdef CFG_CHECK_REPROJECTION
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    if (m_UcsLF[iLF] & LM_FLAG_FRAME_UPDATE_CAMERA_LF) {
      ComputeErrorFeature(&m_LFs[iLF], m_CsLF[iLF].m_T, m_CsKF, m_ds, &m_esLF[iLF].second);
    }
  }
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (m_UcsKF[iKF] & LM_FLAG_FRAME_UPDATE_CAMERA_KF) {
      ComputeErrorFeature(&m_KFs[iKF], m_CsKF[iKF], m_CsKF, m_ds, &m_esKF[iKF].second, iKF);
    }
  }
#endif
  m_LM->LBA_Update(iFrm1, iFrm2, m_ic2LF, m_CsLF, m_UcsLF, m_iFrmsKF, m_CsKF, m_UcsKF,
                   m_iKF2d, m_ds, m_Uds
#ifdef CFG_CHECK_REPROJECTION
                 , m_esLF, m_esKF
#endif
                 );
  const int iLF = m_ic2LF.back(), iFrm = m_LFs[iLF].m_T.m_iFrm;
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_LBA_UpdateData);
  m_C.m_T = m_LFs[iLF].m_T;
  m_C.m_C = m_CsLF[iLF];
  //m_ILFs1.resize(0);
  //m_IKFs1.resize(0);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_LBA_UpdateData);
}

bool LocalBundleAdjustor::BufferDataEmpty() {
  return m_ITs1.empty();
}



void LocalBundleAdjustor::MarginalizeLocalFrame() {
  if (static_cast<int>(m_LFs.size()) < LBA_MAX_LOCAL_FRAMES) {
    return;
  }
  const float eps = 0.0f;
  const float epsr = UT::Inverse(BA_VARIANCE_MAX_ROTATION, BA_WEIGHT_FEATURE, eps);
  const float epsp = UT::Inverse(BA_VARIANCE_MAX_POSITION, BA_WEIGHT_FEATURE, eps);
  const float epsv = UT::Inverse(BA_VARIANCE_MAX_VELOCITY, BA_WEIGHT_FEATURE, eps);
  const float epsba = UT::Inverse(BA_VARIANCE_MAX_BIAS_ACCELERATION, BA_WEIGHT_FEATURE, eps);
  const float epsbw = UT::Inverse(BA_VARIANCE_MAX_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE, eps);
  const float _eps[] = {epsp, epsp, epsp, epsr, epsr, epsr, epsv, epsv, epsv,
                        epsba, epsba, epsba, epsbw, epsbw, epsbw};
  const int iLF1 = m_ic2LF[0], iLF2 = m_ic2LF[1];
  const LocalFrame &LF1 = m_LFs[iLF1];
  const int iKFr = LF1.m_iKFNearest == -1 ? m_Zp.m_iKFr : LF1.m_iKFNearest;

  const int iFrm = m_LFs[m_ic2LF.back()].m_T.m_iFrm;

  m_ZpKF.Invalidate();
  if (LF1.m_T.m_iFrm == m_KFs[iKFr].m_T.m_iFrm) {
    const bool v = m_Zp.Pose::Valid();
    if (v) {
      //m_ts[TM_TEST_2].Start();
      if (m_Zp.GetPriorMotion(&m_ZpLF, &m_work, _eps) ||
          !LBA_MARGINALIZATION_CHECK_INVERTIBLE) {
        //m_GBA->PushCameraPriorMotion(iFrm, m_Zp.m_iKFr, m_ZpLF);
        m_GBA->PushCameraPriorMotion(iFrm, iKFr, m_ZpLF);
      } else {
        //m_ZpLF.MakeZero();
        m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                          BA_VARIANCE_PRIOR_VELOCITY_RESET,
                          BA_VARIANCE_PRIOR_BIAS_ACCELERATION_RESET,
                          BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_RESET, &m_CsLF[iLF1]);
      }
      //m_ts[TM_TEST_2].Stop();
      if ((m_Zp.GetPriorPose(iKFr, &m_ZpKF, &m_work, _eps) ||
           !LBA_MARGINALIZATION_CHECK_INVERTIBLE) &&
           m_ZpKF.MarginalizeUninformative(BA_WEIGHT_FEATURE,
                                           BA_VARIANCE_PRIOR_POSITION_INFORMATIVE,
                                           BA_VARIANCE_PRIOR_ROTATION_INFORMATIVE,
                                           &m_idxsTmp1, &m_work, _eps)) {
        m_GBA->PushCameraPriorPose(iFrm, m_ZpKF);
      } else {
        m_ZpKF.Invalidate();
      }
    }
    const Rigid3D &Tr = m_CsKF[iKFr];
    //const Rigid3D &Tr = m_CsLF[iLF1].m_T;
    const float s2g = v ? BA_VARIANCE_PRIOR_GRAVITY_NEW : BA_VARIANCE_PRIOR_GRAVITY_FIRST;
    //m_ZpLF.m_Amm *= sp;
    //m_ZpLF.m_bm *= sp;
    m_Zp.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL/* * sp*/, iKFr, Tr, s2g, m_ZpLF, true);

  
    Camera C1;
    C1.m_T = Tr;
    Tr.GetPosition(C1.m_p);
    m_Zp.GetMotion(Tr, &C1.m_v, &C1.m_ba, &C1.m_bw);

    IMU::Delta::Error e;
    IMU::Delta::Jacobian::RelativeKF J;
    IMU::Delta::Factor::Auxiliary::RelativeKF A;
    const Camera &C2 = m_CsLF[iLF2];
    const IMU::Delta &D = m_DsLF[iLF2];
    D.GetFactor(BA_WEIGHT_IMU/* * sd*/, C1, C2, m_K.m_pu, &e, &J, &A, BA_ANGLE_EPSILON);

    if (!m_Zp.PropagateKF(Tr, C2, A, &m_work, _eps) && LBA_MARGINALIZATION_CHECK_INVERTIBLE) {
      m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                        BA_VARIANCE_PRIOR_VELOCITY_RESET,
                        BA_VARIANCE_PRIOR_BIAS_ACCELERATION_RESET,
                        BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_RESET, &C2);

      m_Zp.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL, iKFr, Tr,
                      BA_VARIANCE_PRIOR_GRAVITY_RESET, m_ZpLF, false, &C2.m_T,
                      BA_VARIANCE_PRIOR_POSITION_RESET, BA_VARIANCE_PRIOR_ROTATION_RESET);

    }
  } else {
    //Camera C1;
    Rigid3D /*Tr, TrI, */Tr1, Trk, Tk1[2];
    /*const */int Nk = static_cast<int>(m_Zp.m_iKFs.size()) - 1;

    //m_Zp.GetReferencePose(m_CsKF[m_Zp.m_iKFr], &Tr, &TrI);
    //m_Zp.GetPose(TrI, Nk, &Tr1, &C1.m_T);
    m_Zp.m_Zps[Nk].GetInverse(Tr1);
    const float s2d = UT::Inverse(BA_VARIANCE_PRIOR_DEPTH_NEW, BA_WEIGHT_PRIOR_CAMERA_INITIAL); 
    const float epsd = UT::Inverse(BA_VARIANCE_MAX_DEPTH, BA_WEIGHT_FEATURE, eps);
    const float epsc[12] = {epsp, epsp, epsp, epsr, epsr, epsr,
                            epsp, epsp, epsp, epsr, epsr, epsr};
    //const float s2pMax = BA_VARIANCE_PRIOR_POSITION_INSERT_MAX / BA_WEIGHT_FEATURE;
    //const float s2rMax = BA_VARIANCE_PRIOR_ROTATION_INSERT_MAX / BA_WEIGHT_FEATURE;
    //const float r2Max = ME_VARIANCE_HUBER;
    const float r2Max = FLT_MAX;
    const int NZ = static_cast<int>(LF1.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF1.m_Zs[iZ];
      const int iKF = Z.m_iKF;
      const Depth::InverseGaussian *ds = m_ds.data() + m_iKF2d[iKF];
      //const Depth::InverseGaussian *ds = m_dsGT->data() + m_solver->m_internal->m_iKF2d[iKF];
      const KeyFrame &KF = m_KFs[iKF];
      if (iKF == m_Zp.m_iKFr) {
        *Tk1 = Tr1;

        FTR::Factor::FixSource::L L;
        FTR::Factor::FixSource::A1 A1;
        FTR::Factor::FixSource::A2 A2;
        FTR::Factor::FixSource::U U;
        FTR::Factor::FixSource::M2 M;
        FTR::Factor::DD mdd;
        LA::AlignedVector6f mdcz;
        Camera::Factor::Unitary::CC SAczz;
        SAczz.MakeZero();
        for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
          const FTR::Measurement &z = LF1.m_zs[iz];
          const int ix = z.m_ix;
          FTR::GetFactor<LBA_ME_FUNCTION>(BA_WEIGHT_FEATURE, Tk1, KF.m_xs[ix], ds[ix], Tr1, z,
                                          &L, &A1, &A2, &U,
                                          r2Max);
          SAczz += A2.m_Aczz;
          A2.m_add.m_a += s2d;
          if (A2.m_add.m_a < epsd) {
            continue;
          }
          mdd.m_a = 1.0f / A2.m_add.m_a;
          mdd.m_b = mdd.m_a * A2.m_add.m_b;
          FTR::Factor::FixSource::Marginalize(mdd, A1.m_adczA, &mdcz, &M);
          SAczz -= M.m_Mczz;
        }
        m_Zp.Update(Nk, SAczz);
      } else {
        const std::vector<int>::const_iterator it = std::lower_bound(m_Zp.m_iKFs.begin(),
                                                                     m_Zp.m_iKFs.end(), iKF);
        const int ik = static_cast<int>(it - m_Zp.m_iKFs.begin());
        if (it == m_Zp.m_iKFs.end() || *it != iKF) {
          //Rigid3D::ABI(m_CsKF[iKF], Tr, Trk);
          Rigid3D::ABI(m_CsKF[iKF], m_CsKF[m_Zp.m_iKFr], Trk);
        } else {
          m_Zp.m_Zps[ik].GetInverse(Trk);
        }
        Rigid3D::ABI(Tr1, Trk, *Tk1);
        xp128f mdd;
        FTR::Factor::Full::L L;
        FTR::Factor::Full::A1 A1;
        FTR::Factor::Full::A2 A2;
        FTR::Factor::Full::U U;
        FTR::Factor::Full::Source::M1 M1;
        FTR::Factor::Full::Source::M2 M2;
        FTR::Factor::Full::M1 M3;
        FTR::Factor::Full::M2 M4;
        LA::ProductVector6f adcz;
        Camera::Factor::Unitary::CC SAcxx, SAczz;
        Camera::Factor::Binary::CC SAcxz;
        SAcxx.MakeZero();
        SAcxz.MakeZero();
        SAczz.MakeZero();
        for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
          const FTR::Measurement &z = LF1.m_zs[iz];
          const int ix = z.m_ix;
          FTR::GetFactor<LBA_ME_FUNCTION>(BA_WEIGHT_FEATURE, Tk1, KF.m_xs[ix], ds[ix], Tr1, z,
                                          &L, &A1, &A2, &U,
                                          r2Max);
          SAcxx += A2.m_Acxx;
          SAcxz += A2.m_Acxz;
          SAczz += A2.m_Aczz;
          A2.m_adx.m_add.m_a += s2d;

          if (A2.m_adx.m_add.m_a < epsd) {
            continue;
          }
          mdd.vdup_all_lane(1.0f / A2.m_adx.m_add.m_a);
          FTR::Factor::Full::Source::Marginalize(mdd, A2.m_adx, &M1, &M2);
          FTR::Factor::Full::Marginalize(mdd, M1, A1, &M3, &M4, &adcz);
          SAcxx -= M2.m_Mcxx;
          SAcxz -= M4.m_Mcxz;
          SAczz -= M4.m_Mczz;
        }
        if (it == m_Zp.m_iKFs.end() || *it != iKF) {
          if (LBA_MARGINALIZATION_CHECK_RANK) {

            LA::AlignedMatrix12x12f A;
            A.Set(SAcxx.m_A, SAcxz, SAczz.m_A);
            //if (!A.DecomposeLDL(epsc)) {
            //  continue;
            //}
            const int rank = A.RankLDL(epsc);
            //UT::Print("[%d] --> [%d] %d\n", m_KFs[Z.m_iKF].m_T.m_iFrm, LF1.m_T.m_iFrm, rank);
            if (rank < 6) {
            //if (rank < 12) {
              continue;
            }

          }

          const CameraPrior::Element::CC AczzBkp = m_Zp.m_Acc[Nk][Nk];
          const CameraPrior::Element::C bczBkp = m_Zp.m_bc[Nk];
          m_Zp.Insert(BA_WEIGHT_PRIOR_CAMERA_INITIAL, ik, iKF, Trk,
                      BA_VARIANCE_PRIOR_POSITION_NEW,
                      BA_VARIANCE_PRIOR_ROTATION_NEW, &m_work);
          ++Nk;
          m_Zp.Update(ik, Nk, SAcxx, SAcxz, SAczz);

        } else {
          m_Zp.Update(ik, Nk, SAcxx, SAcxz, SAczz);
        }

      }
    }

    if (iKFr == m_Zp.m_iKFr) {
      //C1.m_T.GetPosition(C1.m_p);
    } else {
      //Tr = m_CsKF[iKFr];
      //C1 = m_CsLF[iLF1];
      const Camera &C1 = m_CsLF[iLF1];
      //m_ts[TM_TEST_2].Start();
      if (!m_Zp.GetPriorMotion(&m_ZpLF, &m_work, _eps) && LBA_MARGINALIZATION_CHECK_INVERTIBLE) {
        m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                          BA_VARIANCE_PRIOR_VELOCITY_RESET,
                          BA_VARIANCE_PRIOR_BIAS_ACCELERATION_RESET,
                          BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_RESET, &C1);
#ifdef LBA_DEBUG_GROUND_TRUTH_MEASUREMENT
        if (m_CsGT) {
          m_ZpLF.DebugSetMeasurement(m_CsLFGT[iLF1]);
        }
#endif
      }
      //m_ts[TM_TEST_2].Stop();
      if ((m_Zp.GetPriorPose(INT_MAX, &m_ZpKF, &m_work, _eps) ||
           !LBA_MARGINALIZATION_CHECK_INVERTIBLE) &&
           m_ZpKF.MarginalizeUninformative(BA_WEIGHT_FEATURE,
                                           BA_VARIANCE_PRIOR_POSITION_INFORMATIVE,
                                           BA_VARIANCE_PRIOR_ROTATION_INFORMATIVE,
                                           &m_idxsTmp1, &m_work, _eps)) {
        m_GBA->PushCameraPriorPose(iFrm, m_ZpKF);
      } else {
        m_ZpKF.Invalidate();
      }
      const Rigid3D &Tr = m_CsKF[iKFr];

      m_Zp.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL, iKFr, Tr,
                      BA_VARIANCE_PRIOR_GRAVITY_NEW, m_ZpLF, false, &C1.m_T,
                      BA_VARIANCE_PRIOR_POSITION_NEW, BA_VARIANCE_PRIOR_ROTATION_NEW);
      Rigid3D::ABI(C1.m_T, Tr, Tr1);

    }
    //const Camera &C2 = m_CsLF[iLF2];
    const Rigid3D &Tr = m_CsKF[iKFr];
    Rigid3D _Tr, TrI, Tr2;
    Camera C1, C2;
    LA::AlignedVector3f v2;
    m_Zp.GetReferencePose(Tr, &_Tr, &TrI);
    Rigid3D::ABI(Tr1, TrI, C1.m_T);
    C1.m_T.GetPosition(C1.m_p);
    m_Zp.GetMotion(C1.m_T, &C1.m_v, &C1.m_ba, &C1.m_bw);
    C2 = m_CsLF[iLF2];
    Rigid3D::ABI(C2.m_T, Tr, Tr2);
    C2.m_T.ApplyRotation(C2.m_v, v2);
    Rigid3D::ABI(Tr2, TrI, C2.m_T);
    C2.m_T.GetPosition(C2.m_p);
    C2.m_T.ApplyRotationInversely(v2, C2.m_v);

    IMU::Delta::Error e;
    IMU::Delta::Jacobian::RelativeLF J;
    IMU::Delta::Factor::Auxiliary::RelativeLF A;
    const IMU::Delta &D = m_DsLF[iLF2];
    D.GetFactor(BA_WEIGHT_IMU/* * sd*/, C1, C2, m_K.m_pu, _Tr, &e, &J, &A, BA_ANGLE_EPSILON);


    if (!m_Zp.PropagateLF(_Tr, C2, A, &m_work, _eps) && LBA_MARGINALIZATION_CHECK_INVERTIBLE) {
      const Camera &_C2 = m_CsLF[iLF2];
      m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                        BA_VARIANCE_PRIOR_VELOCITY_RESET,
                        BA_VARIANCE_PRIOR_BIAS_ACCELERATION_RESET,
                        BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_RESET, &_C2);

      m_Zp.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL, iKFr, Tr,
                      BA_VARIANCE_PRIOR_GRAVITY_RESET, m_ZpLF, false, &_C2.m_T,
                      BA_VARIANCE_PRIOR_POSITION_RESET, BA_VARIANCE_PRIOR_ROTATION_RESET);

    }

  }
//#if 0
#if 1
  const int Nk = static_cast<int>(m_Zp.m_iKFs.size() - 1);

  m_Zp.m_Acc[Nk][Nk].IncreaseDiagonal(UT::Inverse(BA_VARIANCE_PRIOR_POSITION_NEW,
                                                  BA_WEIGHT_PRIOR_CAMERA_INITIAL),
                                      UT::Inverse(BA_VARIANCE_PRIOR_ROTATION_NEW,
                                                  BA_WEIGHT_PRIOR_CAMERA_INITIAL));
  m_Zp.m_Amm.IncreaseDiagonal(UT::Inverse(BA_VARIANCE_PRIOR_VELOCITY_NEW,
                                          BA_WEIGHT_PRIOR_CAMERA_INITIAL),
                              UT::Inverse(BA_VARIANCE_PRIOR_BIAS_ACCELERATION_NEW,
                                          BA_WEIGHT_PRIOR_CAMERA_INITIAL),
                              UT::Inverse(BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_NEW,
                                          BA_WEIGHT_PRIOR_CAMERA_INITIAL));
#endif
  //m_ts[TM_TEST_2].Start();
  if (!m_Zp.GetPriorMotion(&m_ZpLF, &m_work, _eps) && LBA_MARGINALIZATION_CHECK_INVERTIBLE) {
    m_ZpLF.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                      BA_VARIANCE_PRIOR_VELOCITY_RESET,
                      BA_VARIANCE_PRIOR_BIAS_ACCELERATION_RESET,
                      BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_RESET, &m_CsLF[iLF2]);

  }
  //m_ts[TM_TEST_2].Stop();
#ifdef CFG_CAMERA_PRIOR_SQUARE_FORM
  CameraPrior::Element::MM Amm = m_ZpLF.m_Amm;
  CameraPrior::Element::M bm = m_ZpLF.m_bm;
  if (Amm.SolveLDL(bm, _eps + 6)) {
    m_ZpLF.m_em.Set(bm);
  } else {
    m_ZpLF.m_Amm.MakeZero();
    m_ZpLF.m_em.MakeZero();
  }
#endif
  m_ApLF.MakeZero();
  m_ucmsLF[iLF2] |= LBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                    LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                    LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
}

void LocalBundleAdjustor::PopLocalFrame() {

  //m_ts[TM_TEST_1].Start();
  MarginalizeLocalFrame();
  //m_ts[TM_TEST_1].Stop();

  const int nLFs = static_cast<int>(m_LFs.size());
  const int STL = std::min(nLFs, LBA_MAX_SLIDING_TRACK_LENGTH);
  const int iLF = m_ic2LF.front();
  LocalFrame &LF = m_LFs[iLF];
  const int NZ = static_cast<int>(LF.m_Zs.size());
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = LF.m_Zs[iZ];
    m_idxsListTmp.resize(STL);
    for (int i = 0; i < STL; ++i) {
      m_idxsListTmp[i].resize(0);
    }
    bool popST = false;
    const int iKF = Z.m_iKF;
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    const ubyte pushFrm = m_ucsKF[iKF] & LBA_FLAG_FRAME_PUSH_TRACK;
    KeyFrame &KF = m_KFs[iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      const int ix = LF.m_zs[iz].m_ix, iSTMin = KF.m_ix2ST[ix];
      int icMin = -1;
      for (int _ic = 1; _ic < STL && icMin == -1; ++_ic) {
        std::vector<int> &_iz2x = m_idxsListTmp[_ic];
        if (_iz2x.empty()) {
          MarkFeatureMeasurements(m_LFs[m_ic2LF[_ic]], iKF, _iz2x);
        }
        if (_iz2x[ix] != -1) {
          icMin = _ic;
        }
      }
      FTR::Factor::FixSource::A2 &Az = LF.m_Azs2[iz];
      Az.m_add.MakeMinus();
      FTR::Factor::FixSource::A3 &AzST = LF.m_AzsST[iz];
      AzST.m_add.MakeMinus();

      if (icMin == -1 ||
          (KF.m_ix2ST[ix + 1] - iSTMin > 1 && KF.m_STs[iSTMin + 1].m_icMin == icMin)) {
        popST = true;
        uds[ix] |= LBA_FLAG_TRACK_POP;
        KF.m_usST[iSTMin] |= LBA_FLAG_TRACK_POP;
        FTR::Factor::DD &mddST = KF.m_MxsST[iSTMin].m_mdd;
        mddST.MakeMinus();
        const bool nonZero = !(KF.m_usST[iSTMin] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO);
        if (icMin != -1) {
          for (int _ic = 1; _ic < STL; ++_ic) {
            const int _iLF = m_ic2LF[_ic];
            LocalFrame &_LF = m_LFs[_iLF];
            std::vector<int> &_ix2z = m_idxsListTmp[_ic];
            if (_ix2z.empty()) {
              MarkFeatureMeasurements(_LF, iKF, _ix2z);
            }
            const int _iz = _ix2z[ix];
            if (_iz == -1) {
              continue;
            }
            --_LF.m_STs[_iz].m_ist2;
            if (nonZero) {
              _LF.m_SmddsST[_iz] += mddST;
              _LF.m_ms[_iz] |= LBA_FLAG_MARGINALIZATION_UPDATE;
            }
          }
        }
        if (nonZero) {
          KF.m_ms[ix] |= LBA_FLAG_MARGINALIZATION_UPDATE;
        }
      } else {
        KF.m_STs[iSTMin].m_icMin = icMin;
        if (!pushFrm || !(KF.m_usST[iSTMin] & LBA_FLAG_TRACK_PUSH)) {
          KF.m_usST[iSTMin] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
          FTR::Factor::DD &SaddST = KF.m_AxsST[iSTMin].m_Sadd;

          SaddST += AzST.m_add;
          if (SaddST.m_a < 0.0f) {
            //SaddST.MakeZero();
            uds[ix] |= LBA_FLAG_TRACK_UPDATE_DEPTH;
            m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
          }
        }
      }
      FTR::Factor::DD &Sadd = KF.m_Axs[ix].m_Sadd;
      Sadd += Az.m_add;
      if (Sadd.m_a < 0.0f) {
        //Sadd.MakeZero();
        uds[ix] |= LBA_FLAG_TRACK_UPDATE_DEPTH;
        m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
      }
      //m_F = -Az.m_F + m_F;
      uds[ix] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
    }
    if (Z.m_iz1 < Z.m_iz2) {
      m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
    }
    if (popST) {
      m_ucsKF[iKF] |= LBA_FLAG_FRAME_POP_TRACK;
    }
  }
  const int icMax = STL - 1;
  for (int _ic = 1; _ic < icMax; ++_ic) {
    const int iLF1 = m_ic2LF[_ic];
    LocalFrame &LF1 = m_LFs[iLF1];
    const int NI = int(LF1.m_Zm.m_Is.size());
    for (int iI = 0; iI < NI; ++iI) {
      const MeasurementMatchLF::Index &I = LF1.m_Zm.m_Is[iI];
      if (!(m_ucsKF[I.m_iKF] & LBA_FLAG_FRAME_POP_TRACK)) {
        continue;
      }
      const KeyFrame &KF = m_KFs[I.m_iKF];
      const int iLF2 = LF1.m_iLFsMatch[I.m_ik];
      const LocalFrame &LF2 = m_LFs[iLF2];
      const int i1 = LF1.m_Zm.m_iI2zm[iI], i2 = LF1.m_Zm.m_iI2zm[iI + 1];
      for (int i = i1; i < i2; ++i) {
        const int iz2 = LF1.m_Zm.m_izms[i].m_iz2, ix = LF2.m_zs[iz2].m_ix, iST = KF.m_ix2ST[ix];
        if (LF2.m_STs[iz2].m_ist1 != 0 || !(KF.m_usST[iST] & LBA_FLAG_TRACK_POP) ||
            (KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
          continue;
        }
        //LF1.m_Zm.m_SmddsST[i] = -KF.m_MxsST[iST].m_mdd.m_a + LF1.m_ZmLF.m_SmddsST[i];
        LF1.m_Zm.m_SmddsST[i] = KF.m_MxsST[iST].m_mdd.m_a + LF1.m_Zm.m_SmddsST[i];
        LF1.m_Zm.m_ms[i] |= LBA_FLAG_MARGINALIZATION_UPDATE;
      }
    }
  }
  for (int _ic = STL; _ic < nLFs; ++_ic) {
    LocalFrame &_LF = m_LFs[m_ic2LF[_ic]];
    const int _NZ = int(_LF.m_Zs.size());
    for (int _iZ = 0; _iZ < _NZ; ++_iZ) {
      const FRM::Measurement &_Z = _LF.m_Zs[_iZ];
      if (!(m_ucsKF[_Z.m_iKF] & LBA_FLAG_FRAME_POP_TRACK)) {
        continue;
      }
      const KeyFrame &KF = m_KFs[_Z.m_iKF];
      const int _iz1 = _Z.m_iz1, _iz2 = _Z.m_iz2;
      for (int _iz = _iz1; _iz < _iz2; ++_iz) {
        if (KF.m_usST[KF.m_ix2ST[_LF.m_zs[_iz].m_ix]] & LBA_FLAG_TRACK_POP) {
          _LF.m_STs[_iz].Step();
        }
      }
    }
  }
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = LF.m_Zs[iZ];
    if (!(m_ucsKF[Z.m_iKF] & LBA_FLAG_FRAME_POP_TRACK)) {
      continue;
    }
    KeyFrame &KF = m_KFs[Z.m_iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    m_ix2STTmp.swap(KF.m_ix2ST);
    KF.m_ix2ST.resize(Nx + 1);
    KF.m_ix2ST[0] = 0;
    for (int ix = 0, iST = 0; ix < Nx; ++ix) {
      const int iST1 = m_ix2STTmp[ix], iST2 = m_ix2STTmp[ix + 1];
      const bool popST = iST1 < iST2 && (KF.m_usST[iST1] & LBA_FLAG_TRACK_POP);
      for (int _iST = popST ? iST1 + 1 : iST1; _iST < iST2; ++iST, ++_iST) {
        KF.m_STs[iST] = KF.m_STs[_iST];
        KF.m_usST[iST] = KF.m_usST[_iST];
        KF.m_AxsST[iST] = KF.m_AxsST[_iST];
        KF.m_MxsST[iST] = KF.m_MxsST[_iST];
      }
      KF.m_ix2ST[ix + 1] = iST;
      if (popST && iST == KF.m_ix2ST[ix]) {
        KF.m_Axs[ix] = KF.m_Axps[ix];
        KF.m_ms[ix] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
      }
    }
    const int NST = KF.m_ix2ST.back();
    if (NST == 0) {
      KF.m_STs.clear();
      KF.m_usST.clear();
      KF.m_AxsST.Clear();
      KF.m_MxsST.Clear();
    } else {
      KF.m_STs.resize(NST);
      KF.m_usST.resize(NST);
      KF.m_AxsST.Resize(NST);
      KF.m_MxsST.Resize(NST);
    }
  }
  const int _iLF = m_ic2LF[1];
  IMU::Delta::Factor &Ad = m_AdsLF[_iLF];
  Ad.m_A22.MakeMinus();
  m_SAcusLF[_iLF] += Ad.m_A22.m_Acc;
  Camera::Factor &SAcm = m_SAcmsLF[_iLF];
  SAcm.m_Au.m_Acm += Ad.m_A22.m_Acm;
  SAcm.m_Au.m_Amm += Ad.m_A22.m_Amm;
  SAcm.m_Ab.MakeZero();
  Ad.MakeZero();
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    KeyFrame &KF = m_KFs[iKF];
    const int NST = static_cast<int>(KF.m_STs.size());
    for (int iST = 0; iST < NST; ++iST) {
      KF.m_STs[iST].Step();
    }
  }
}

void LocalBundleAdjustor::_PushLocalFrame(const InputLocalFrame &ILF) {
  const int nLFs1 = static_cast<int>(m_LFs.size());
  if (static_cast<int>(m_LFs.size()) < LBA_MAX_LOCAL_FRAMES) 
  {
    const int nLFs2 = nLFs1 + 1;
    m_ic2LF.push_back(nLFs1);
    m_LFs.resize(nLFs2);
    m_CsLF.Resize(nLFs2, true);

    m_ucsLF.resize(nLFs2);
    m_ucmsLF.resize(nLFs2);
#ifdef CFG_INCREMENTAL_PCG
    m_xcsLF.Resize(nLFs2, true);
    m_xmsLF.Resize(nLFs2, true);
#endif
    m_DsLF.Resize(nLFs2, true);

    m_AdsLF.Resize(nLFs2, true);
    m_AfpsLF.Resize(nLFs2, true);
    m_AfmsLF.Resize(nLFs2, true);
    m_SAcusLF.Resize(nLFs2, true);
    m_SMcusLF.Resize(nLFs2, true);
    m_SAcmsLF.Resize(nLFs2, true);
    m_UcsLF.resize(nLFs2, LM_FLAG_FRAME_DEFAULT);
  } else {
    PopLocalFrame();
    const int iLF = m_ic2LF.front();
    m_ic2LF.erase(m_ic2LF.begin());
    m_ic2LF.push_back(iLF);
  }

  const int iLF = m_ic2LF.back();
  LocalFrame &LF = m_LFs[iLF];

  LF.Initialize(ILF, ILF.m_us);
  while (!LF.m_Zs.empty() && LF.m_Zs.back().m_iKF >= static_cast<int>(m_KFs.size())) {
    LF.PopFrameMeasurement();
  }

  std::vector<int> &iKF2X = m_idxsTmp1, &iX2z = m_idxsTmp2;
  PushFeatureMeasurementMatchesFirst(LF, iKF2X, iX2z);
  const int nLFs2 = int(m_LFs.size()), STL = std::min(nLFs2, LBA_MAX_SLIDING_TRACK_LENGTH);
  const int ic1 = nLFs2 - STL, ic2 = nLFs2 - 1;
  for (int _ic = ic1; _ic < ic2; ++_ic) {
    LocalFrame &_LF = m_LFs[m_ic2LF[_ic]];
    _LF.m_iLFsMatch.push_back(iLF);
    PushFeatureMeasurementMatchesNext(_LF, LF, iKF2X, iX2z, _LF.m_Zm);
  }
  m_ucsLF[iLF] = LBA_FLAG_FRAME_UPDATE_CAMERA;
  m_ucmsLF[iLF] = LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION | LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION |
                  LBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                  LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION | LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
  m_UcsLF[iLF] = LM_FLAG_FRAME_UPDATE_CAMERA_LF;

  const ubyte udFlagST = LBA_FLAG_TRACK_PUSH | LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
  const int NZ = static_cast<int>(LF.m_Zs.size());

  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = LF.m_Zs[iZ];
    bool pushST = false;
    ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
    KeyFrame &KF = m_KFs[Z.m_iKF];
    m_marksTmp1.assign(KF.m_xs.size(), 0);

    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      const int ix = LF.m_zs[iz].m_ix, Nst = KF.CountSlidingTracks(ix), iSTMax = KF.m_ix2ST[ix + 1] - 1;
      
      if (Nst > 0 && KF.m_STs[iSTMax].m_icMin >= ic1) {
        LF.m_STs[iz].Set(Nst - 1, Nst);
        KF.m_STs[iSTMax].m_icMax = ic2;
      } else {
        LF.m_STs[iz].Set(Nst, Nst + 1);
        pushST = true;
        uds[ix] |= LBA_FLAG_TRACK_PUSH;
        //uds[ix] |= LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
        m_marksTmp1[ix] = 1;
      }
    }
    if (!pushST) {
      continue;
    }
    m_ucsKF[Z.m_iKF] |= LBA_FLAG_FRAME_PUSH_TRACK;
    m_idxsListTmp.resize(STL);

    for (int i = 0; i < STL; ++i) {
      m_idxsListTmp[i].resize(0);
    }
    const int Nx = static_cast<int>(KF.m_xs.size());
    m_ix2STTmp.swap(KF.m_ix2ST);  KF.m_ix2ST.resize(Nx + 1);
    m_STsTmp.swap(KF.m_STs);      KF.m_STs.resize(0);
    m_usSTTmp.swap(KF.m_usST);    KF.m_usST.resize(0);
    m_AxsTmp.Swap(KF.m_AxsST);    KF.m_AxsST.Resize(0);
    m_MxsTmp.Swap(KF.m_MxsST);    KF.m_MxsST.Resize(0);

    for (int ix = 0; ix < Nx; ++ix) {
      const int iST1 = m_ix2STTmp[ix], iST2 = m_ix2STTmp[ix + 1], Nst = iST2 - iST1;
      KF.m_ix2ST[ix] = static_cast<int>(KF.m_STs.size());
      KF.m_STs.insert(KF.m_STs.end(), m_STsTmp.begin() + iST1, m_STsTmp.begin() + iST2);
      KF.m_usST.insert(KF.m_usST.end(), m_usSTTmp.begin() + iST1, m_usSTTmp.begin() + iST2);
      KF.m_AxsST.Push(m_AxsTmp.Data() + iST1, Nst);
      KF.m_MxsST.Push(m_MxsTmp.Data() + iST1, Nst);
      //if (!(uds[ix] & LBA_FLAG_TRACK_PUSH))
      if (!m_marksTmp1[ix]) {
        continue;
      }
      const int NST1 = static_cast<int>(KF.m_STs.size()), NST2 = NST1 + 1;
      KF.m_STs.resize(NST2);
      KF.m_usST.push_back(udFlagST);
      KF.m_AxsST.Resize(NST2, true);  KF.m_AxsST[NST1].MakeZero();
      KF.m_MxsST.Resize(NST2, true);  KF.m_MxsST[NST1].MakeZero();

      KeyFrame::SlidingTrack &ST = KF.m_STs[NST1];
      ST.Set(ic2);
      for (int _ic = ic1; _ic < ic2; ++_ic) {
        const int _iLF = m_ic2LF[_ic];
        LocalFrame &_LF = m_LFs[_iLF];
        std::vector<int> &_ix2z = m_idxsListTmp[_ic - ic1];
        if (_ix2z.empty()) {
          MarkFeatureMeasurements(_LF, Z.m_iKF, _ix2z);
        }
        const int _iz = _ix2z[ix];
        if (_iz == -1) {
          continue;
        }
        ST.m_icMin = std::min(ST.m_icMin, _ic);

        ++_LF.m_STs[_iz].m_ist2;
      }
    }
    KF.m_ix2ST[Nx] = static_cast<int>(KF.m_STs.size());
  }
  m_CsLF[iLF] = ILF.m_C;

#ifdef CFG_INCREMENTAL_PCG
  m_xcsLF[iLF].MakeZero();
  m_xmsLF[iLF].MakeZero();
#endif
  IMU::Delta &D = m_DsLF[iLF];
  if (nLFs2 > 1) {
    const int _iLF = m_ic2LF[nLFs2 - 2];
    const LocalFrame &_LF = m_LFs[_iLF];
    const float _t = _LF.m_T.m_t;
    //UT::DebugStart();
    IMU::PreIntegrate(LF.m_us, _t, LF.m_T.m_t, m_CsLF[_iLF], &D, &m_work, true,
                      &_LF.m_us.Back(), NULL, BA_ANGLE_EPSILON);
    //UT::DebugStop();
    if (nLFs2 > 2) {
      const int __iLF = m_ic2LF[nLFs2 - 3];
      IMU::Delta &_D = m_DsLF[_iLF];
      IMU::PreIntegrate(_LF.m_us, m_LFs[__iLF].m_T.m_t, _t, m_CsLF[__iLF], &_D, &m_work, true,
                        &_D.m_u1, &LF.m_us.Front(), BA_ANGLE_EPSILON);
      m_ucsLF[_iLF] |= LBA_FLAG_FRAME_UPDATE_CAMERA;
      m_ucmsLF[_iLF] |= LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION | LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION |
                        LBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                        LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION | LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
    }
  } else {
    D.Invalidate();
  }

  if (LF.m_T.m_iFrm == 0) {    
    const LA::Vector3f s2r = LA::Vector3f::Get(BA_VARIANCE_FIX_ORIGIN_ROTATION_X,
                                               BA_VARIANCE_FIX_ORIGIN_ROTATION_Y,
                                               BA_VARIANCE_FIX_ORIGIN_ROTATION_Z);
    const float s2p = BA_VARIANCE_FIX_ORIGIN_POSITION;
    m_Zo.Set(BA_WEIGHT_FIX_ORIGIN, s2r, s2p, m_CsLF[iLF].m_T);
    m_Ao.MakeZero();
  }
  m_AdsLF[iLF].MakeZero();
  m_AfpsLF[iLF].MakeZero();
  m_AfmsLF[iLF].MakeZero();
  m_SAcusLF[iLF].MakeZero();
  m_SMcusLF[iLF].MakeZero();
  m_SAcmsLF[iLF].MakeZero();
  m_usKF.Push(ILF.m_us);

}

void LocalBundleAdjustor::_PushKeyFrame(const GlobalMap::InputKeyFrame &IKF) {
  //Timer timer;
  //timer.Start();
  const int nKFs1 = static_cast<int>(m_KFs.size()), nKFs2 = nKFs1 + 1;
  m_KFs.resize(nKFs2);
  m_iFrmsKF.resize(nKFs2);
  m_CsKF.Resize(nKFs2, true);
  m_ucsKF.resize(nKFs2, LBA_FLAG_FRAME_UPDATE_CAMERA);
  m_UcsKF.resize(nKFs2, LM_FLAG_FRAME_UPDATE_CAMERA_KF);

  const int iKF = nKFs1;
  KeyFrame &KF = m_KFs[iKF];
  KF.Initialize(IKF);
  m_iFrmsKF[iKF] = KF.m_T.m_iFrm;
  //KF.m_zs[368].m_z.Print(true);
  const int Nk = static_cast<int>(KF.m_iKFsMatch.size());

  for (int ik = 0; ik < Nk; ++ik) {
    KeyFrame &_KF = m_KFs[KF.m_iKFsMatch[ik]];

    _KF.m_iKFsMatch.push_back(iKF);
  }
  m_iKF2d.push_back(m_iKF2d.back());
  m_dsBkp.resize(KF.m_zs.size());
  const int NZ1 = static_cast<int>(KF.m_Zs.size());
  for (int iZ = 0; iZ < NZ1; ++iZ) {
    const FRM::Measurement &Z = KF.m_Zs[iZ];
    const Depth::InverseGaussian *ds = m_ds.data() + m_iKF2d[Z.m_iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      m_dsBkp[iz] = ds[KF.m_zs[iz].m_ix];
    }
  }
  const ubyte udFlag1 = LBA_FLAG_TRACK_UPDATE_DEPTH | LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
  const int NX = static_cast<int>(IKF.m_Xs.size());
  for (int iX1 = 0, iX2 = 0; iX1 < NX; iX1 = iX2) {
    const int _iKF = IKF.m_Xs[iX1].m_iKF;
    for (iX2 = iX1 + 1; iX2 < NX && IKF.m_Xs[iX2].m_iKF == _iKF; ++iX2) {}
    const int id = m_iKF2d[_iKF + 1], Nx = iX2 - iX1;
    for (int jKF = _iKF; jKF <= iKF; ++jKF) {
      m_iKF2d[jKF + 1] += Nx;
    }
    m_ds.insert(m_ds.begin() + id, Nx, Depth::InverseGaussian());
    m_uds.insert(m_uds.begin() + id, Nx, udFlag1);
    m_Uds.insert(m_Uds.begin() + id, Nx, LM_FLAG_TRACK_UPDATE_DEPTH);
    m_ucsKF[_iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
    m_UcsKF[_iKF] |= LM_FLAG_FRAME_UPDATE_DEPTH;

    const GlobalMap::Point *Xs = IKF.m_Xs.data() + iX1;
    Depth::InverseGaussian *ds = m_ds.data() + id;
    ubyte *uds = m_uds.data() + id;
    m_xsTmp.resize(Nx);
    std::vector<ubyte> &mcs = m_marksTmp1;
    mcs.assign(nKFs2, 0);
    for (int i = 0; i < Nx; ++i) {
      const GlobalMap::Point &X = Xs[i];
      ds[i] = X.m_d;
      m_xsTmp[i] = X.m_x;
      const int Nz = static_cast<int>(X.m_zs.size());
#ifdef LBA_FLAG_TRACK_MEASURE_KF
      if (Nz > 0) {
        uds[i] |= LBA_FLAG_TRACK_MEASURE_KF;
      }
#endif
      for (int j = 0; j < Nz; ++j) {
        mcs[X.m_zs[j].m_iKF] = 1;
      }
    }

    std::vector<int> &ik2KF = m_idxsTmp1, &iKF2k = m_idxsTmp2;
    ik2KF.resize(0);
    iKF2k.assign(nKFs2, -1);
    for (int jKF = 0; jKF < nKFs2; ++jKF) {
      if (!mcs[jKF]) {
        continue;
      }
      iKF2k[jKF] = static_cast<int>(ik2KF.size());
      ik2KF.push_back(jKF);
    }
    const int Nk = static_cast<int>(ik2KF.size());
    m_zsListTmp.resize(Nk);
    for (int ik = 0; ik < Nk; ++ik) {
      m_zsListTmp[ik].resize(0);
    }
    KeyFrame &_KF = m_KFs[_iKF];
    const int Nx1 = static_cast<int>(_KF.m_xs.size());
    for (int i = 0, ix = Nx1; i < Nx; ++i, ++ix) {
      const GlobalMap::Point &X = Xs[i];
      const int Nz = static_cast<int>(X.m_zs.size());
      for (int j = 0; j < Nz; ++j) {
        const FTR::Measurement &z = X.m_zs[j];
        const int ik = iKF2k[z.m_iKF];
        m_zsListTmp[ik].push_back(z);
        m_zsListTmp[ik].back().m_ix = ix;
      }
    }
    std::vector<int> &iks = m_idxsTmp2;
    iks.resize(Nk);
    std::vector<ubyte> &mxs = m_marksTmp2;
    for (int ik2 = 0; ik2 < Nk; ++ik2) {
      const int iKF2 = ik2KF[ik2];
      KeyFrame &KF2 = m_KFs[iKF2];
      int &jk2 = iks[ik2];
      const std::vector<FTR::Measurement> &zs2 = m_zsListTmp[ik2];
      const int Nk2 = static_cast<int>(KF2.m_iKFsMatch.size());
      KF2.PushFeatureMeasurements(_iKF, zs2, &jk2, &m_work);
      if (static_cast<int>(KF2.m_iKFsMatch.size()) > Nk2) {
        _KF.InsertMatchKeyFrame(iKF2);
      }
      mxs.assign(Nx, -1);
      ubyte *_mxs = mxs.data() - Nx1;
      const int Nz2 = static_cast<int>(zs2.size());
      for (int i = 0; i < Nz2; ++i) {
        const int ix = zs2[i].m_ix;

        _mxs[ix] = 1;
      }
      for (int ik1 = 0; ik1 < ik2; ++ik1) {
        const int iKF1 = ik2KF[ik1];
        const std::vector<FTR::Measurement> &zs1 = m_zsListTmp[ik1];
        const int Nz1 = static_cast<int>(zs1.size());
        bool found = false;
        for (int i = 0; i < Nz1; ++i) {
          const int ix = zs1[i].m_ix;

          if (_mxs[ix]) {
            found = true;
            break;
          }
        }
        if (!found) {
          continue;
        }
        const std::vector<int>::iterator _jk2 = std::lower_bound(KF2.m_iKFsMatch.begin() + jk2,
                                                                 KF2.m_iKFsMatch.end(), iKF1);
        jk2 = static_cast<int>(_jk2 - KF2.m_iKFsMatch.begin());
        if (_jk2 != KF2.m_iKFsMatch.end() && *_jk2 == iKF1) {
          continue;
        }
        KF2.InsertMatchKeyFrame(iKF1, &_jk2);
        KeyFrame &KF1 = m_KFs[iKF1];
        int &jk1 = iks[ik1];
        const std::vector<int>::iterator _jk1 = std::lower_bound(KF1.m_iKFsMatch.begin() + jk1,
                                                                 KF1.m_iKFsMatch.end(), iKF2);
        jk1 = static_cast<int>(_jk1 - KF1.m_iKFsMatch.begin());
        KF1.InsertMatchKeyFrame(iKF2, &_jk1);
      }
    }
    _KF.PushFeatures(m_xsTmp);
  }
  std::vector<int> &iKF2X = m_idxsTmp1, &iX2z = m_idxsTmp2;
  PushFeatureMeasurementMatchesFirst(KF, iKF2X, iX2z);
  m_CsKF[iKF] = IKF.m_C.m_T;


  const ubyte udFlag2 = LBA_FLAG_TRACK_PUSH | LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
  const int nLFs = static_cast<int>(m_ic2LF.size());
  for (int ic = 0; ic < nLFs; ++ic) {
    const int iLF = m_ic2LF[ic];
    LocalFrame &LF = m_LFs[iLF];

    if (LF.m_T == KF.m_T) {
      const int Nx = static_cast<int>(KF.m_xs.size()), Nz = static_cast<int>(LF.m_zs.size());

      if (Nx > 0) {
        LF.PushFrameMeasurement(nKFs1, Nx);
      } else {
        bool found = false;
        const int NZ = static_cast<int>(LF.m_Zs.size());
        for (int iZ = 0; iZ < NZ && !found; ++iZ) {
          const FRM::Measurement &Z = LF.m_Zs[iZ];
          const int iX = iKF2X[Z.m_iKF];
          if (iX == -1) {
            continue;
          }
          const int *_ix2z = iX2z.data() + iX;
          for (int iz = Z.m_iz1; iz < Z.m_iz2 && !found; ++iz) {
            found = _ix2z[LF.m_zs[iz].m_ix] != -1;
          }
        }
        if (found) {
          LF.m_iKFsMatch.push_back(nKFs1);
        }
      }
      KF.m_iKFNearest = LF.m_T.m_iFrm == 0 ? 0 : LF.m_iKFNearest;
      LF.m_iKFNearest = nKFs1;
      m_ucsKF[nKFs1] |= LBA_FLAG_FRAME_PUSH_TRACK;
      ubyte *uds = m_uds.data() + m_iKF2d[nKFs1];
      const GlobalMap::Point *Xs = IKF.m_Xs.data() + IKF.m_Xs.size() - Nx;

      for (int ix = 0, iz = Nz; ix < Nx; ++ix, ++iz) {
        const FTR::Source &x = KF.m_xs[ix];
        LF.m_zs[iz].Set(ix, x.m_x, Xs[ix].m_W
                      );
        LF.m_STs[iz].Set(0, 1);
        KF.m_ix2ST[ix] = ix;
        uds[ix] |= udFlag2;
      }
      KF.m_ix2ST[Nx] = Nx;
      KF.m_STs.assign(Nx, KeyFrame::SlidingTrack(ic));
      KF.m_usST.assign(Nx, udFlag2);
      KF.m_AxsST.Resize(Nx);  KF.m_AxsST.MakeZero();
      KF.m_MxsST.Resize(Nx);  KF.m_MxsST.MakeZero();

    } else {
      bool found = false;
      const int NZ = static_cast<int>(LF.m_Zs.size());
      for (int iZ = 0; iZ < NZ && !found; ++iZ) {
        const FRM::Measurement &Z = LF.m_Zs[iZ];
        const int iX = iKF2X[Z.m_iKF];
        if (iX == -1) {
          continue;
        }
        const int *_ix2z = iX2z.data() + iX;
        for (int iz = Z.m_iz1; iz < Z.m_iz2 && !found; ++iz) {
          found = _ix2z[LF.m_zs[iz].m_ix] != -1;
        }
      }
      if (found) {
        LF.m_iKFsMatch.push_back(nKFs1);
      }
    }
  }
  //const ubyte udFlag3 = LBA_FLAG_TRACK_MEASURE_KF;
  const bool ud = LBA_RESET_DEPTH_INFORMATION;
  const ubyte udFlag3 = (ud ? LBA_FLAG_TRACK_UPDATE_DEPTH : LBA_FLAG_TRACK_DEFAULT)
#ifdef LBA_FLAG_TRACK_MEASURE_KF
                       | LBA_FLAG_TRACK_MEASURE_KF
#endif
                       ;
  const int NZ2 = static_cast<int>(KF.m_Zs.size());
  for (int iZ = 0; iZ < NZ2; ++iZ) {
    const FRM::Measurement &Z = KF.m_Zs[iZ];
    if (ud) {
      m_ucsKF[Z.m_iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
    }
    ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      uds[KF.m_zs[iz].m_ix] |= udFlag3;
    }
  }
  m_GM->LBA_PushKeyFrame(GlobalMap::Camera(IKF.m_C.m_T, KF.m_T.m_iFrm, GM_FLAG_FRAME_DEFAULT
  ));
  const int N = static_cast<int>(std::lower_bound(m_usKF.Data(), m_usKF.End(), KF.m_T.m_t) -
                                                  m_usKF.Data());
  m_usKF.Erase(N, m_usKFLast);
  m_GBA->PushKeyFrame(IKF, m_usKFLast, m_dsBkp

  );

#ifdef CFG_CHECK_REPROJECTION
  m_esKF.resize(nKFs2);
  ComputeErrorFeature(&KF, m_CsKF[nKFs1], m_CsKF, m_ds, &m_esKF[nKFs1].first, nKFs1);
#endif

  //timer.Stop(true);
  //static double g_St = 0.0;
  //static int g_N = 0;
  //const double t = timer.GetAverageMilliseconds();
  //g_St += t;
  //++g_N;
  //UT::Print("[%d] LBA::PushKeyFrame = %f %f\n", IKF.m_T.m_iFrm, t, g_St / g_N);
}

void LocalBundleAdjustor::DeleteKeyFrame(const int iKF) {
  //Timer timer;
  //timer.Start();
  const int iFrm = m_iFrmsKF[iKF];
  const int nKFs1 = static_cast<int>(m_KFs.size()), nKFs2 = nKFs1 - 1;

  const int Nd = static_cast<int>(m_KFs[iKF].m_xs.size());
  const int id1 = m_iKF2d[iKF], id2 = id1 + Nd;
  {
    FTR::Factor::DD daddST;
    KeyFrame &KF = m_KFs[iKF];
    const ubyte ucFlag = LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION |
                         LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION_KF;
    const ubyte udFlag = LBA_FLAG_TRACK_UPDATE_INFORMATION | LBA_FLAG_TRACK_UPDATE_INFORMATION_KF;
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iKF = Z.m_iKF;
      ubyte *_uds = m_uds.data() + m_iKF2d[_iKF];
      KeyFrame &_KF = m_KFs[_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = KF.m_zs[iz].m_ix;
        FTR::Factor::Depth &A = KF.m_Azs[iz];
        A.m_add.MakeMinus();
        _KF.m_Axps[ix].m_Sadd += A.m_add;
        _KF.m_Axs[ix].m_Sadd += A.m_add;
        _uds[ix] |= udFlag;
        FTR::Factor::FixSource::Source::A &AST = _KF.m_AxpsST[ix];
        daddST = AST.m_Sadd;
        AST = _KF.m_Axps[ix];
        const int iST1 = _KF.m_ix2ST[ix], iST2 = _KF.m_ix2ST[ix + 1], Nst = iST2 - iST1;
        if (Nst > 1) {
          _KF.m_AxpsST[ix] *= 1.0f / Nst;
        }
        FTR::Factor::DD::amb(AST.m_Sadd, daddST, daddST);
        for (int iST = iST1; iST < iST2; ++iST) {
          _KF.m_AxsST[iST].m_Sadd += daddST;
          _KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
        }
      }
      if (Z.m_iz1 < Z.m_iz2) {
        m_ucsKF[_iKF] |= ucFlag;
      }
    }
  }
  for (int jKF = 0; jKF < iKF; ++jKF) {
    m_KFs[jKF].DeleteMatchKeyFrame(iKF);
  }
  for (int jKF = iKF + 1; jKF < nKFs1; ++jKF) {
    m_KFs[jKF].DeleteKeyFrame(iKF);
  }
  m_KFs.erase(m_KFs.begin() + iKF);
  m_iFrmsKF.erase(m_iFrmsKF.begin() + iKF);
  for (int jKF = iKF + 1; jKF <= nKFs1; ++jKF) {
    m_iKF2d[jKF] -= Nd;
  }
  m_iKF2d.erase(m_iKF2d.begin() + iKF);
  m_CsKF.Erase(iKF);
  m_ucsKF.erase(m_ucsKF.begin() + iKF);
  m_UcsKF.erase(m_UcsKF.begin() + iKF);

  m_ds.erase(m_ds.begin() + id1, m_ds.begin() + id2);
  m_uds.erase(m_uds.begin() + id1, m_uds.begin() + id2);
  //m_Uds.erase(m_Uds.begin() + id1, m_Uds.begin() + id2);
  const int nLFs = static_cast<int>(m_LFs.size());
  m_iLF2Z.resize(nLFs);
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    LocalFrame &LF = m_LFs[iLF];
    m_iLF2Z[iLF] = std::lower_bound(LF.m_Zs.begin(), LF.m_Zs.end(), iKF);
  }
  for (int iLF1 = 0; iLF1 < nLFs; ++iLF1) {
    LocalFrame &LF1 = m_LFs[iLF1];
    const std::vector<FRM::Measurement>::iterator iZ1 = m_iLF2Z[iLF1];
    const bool z1 = iZ1 != LF1.m_Zs.end() && iZ1->m_iKF == iKF;
    const int Nk = static_cast<int>(LF1.m_iLFsMatch.size());
    for (int ik = 0; ik < Nk; ++ik) {
      const int iLF2 = LF1.m_iLFsMatch[ik];
      const std::vector<FRM::Measurement>::iterator iZ2 = m_iLF2Z[iLF2];
      const bool z2 = iZ2 != m_LFs[iLF2].m_Zs.end() && iZ2->m_iKF == iKF;
      if (!z1 && !z2) {
        continue;
      }
      std::vector<FTR::Measurement::Match> &izms = LF1.m_Zm.m_izms;
      const int i1 = LF1.m_Zm.m_ik2zm[ik], i2 = LF1.m_Zm.m_ik2zm[ik + 1];
      if (z1) {
        const int Nz1 = iZ1->CountFeatureMeasurements();
        for (int i = i2 - 1; i >= i1 && izms[i].m_iz1 >= iZ1->m_iz2; --i) {
          izms[i].m_iz1 -= Nz1;
        }
      }
      if (z2) {
        const int Nz2 = iZ2->CountFeatureMeasurements();
        for (int i = i2 - 1; i >= i1 && izms[i].m_iz2 >= iZ2->m_iz2; --i) {
          izms[i].m_iz2 -= Nz2;
        }
      }
    }
  }
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    LocalFrame &LF = m_LFs[iLF];
    const std::vector<FRM::Measurement>::iterator iZ = m_iLF2Z[iLF];
    if (iZ != LF.m_Zs.end() && iZ->m_iKF == iKF) {
      Camera::Factor::Unitary::CC &SAczz = m_SAcusLF[iLF], &SMczz = m_SMcusLF[iLF];
      for (int iz = iZ->m_iz1; iz < iZ->m_iz2; ++iz) {
        Camera::Factor::Unitary::CC &Aczz = LF.m_Azs2[iz].m_Aczz;
        Aczz.MakeMinus();
        SAczz += Aczz;
        if (!(LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_NON_ZERO)) {
          continue;
        }
        Camera::Factor::Unitary::CC &Mczz = LF.m_Mzs2[iz].m_Mczz;
        Mczz.MakeMinus();
        SMczz += Mczz;
      }
    }
    LF.DeleteKeyFrame(iKF, &iZ);
    if (LF.m_iKFNearest != -1) {
      continue;
    }
    if (LBA_MARGINALIZATION_REFERENCE_NEAREST) {
      ubyte first = 1;
      int iKFNearest = -1;
      float imgMotionNearest = FLT_MAX;
      const Rigid3D &C = m_CsLF[iLF].m_T;
      const float z = 1.0f / LF.m_d.u();
      m_marksTmp1.assign(nKFs2, 0);
      const int Nk = static_cast<int>(LF.m_iKFsMatch.size());
      for (int i = 0; i < Nk; ++i) {
        m_marksTmp1[LF.m_iKFsMatch[i]] = 1;
      }
      for (int jKF = 0; jKF < nKFs2 && m_KFs[jKF].m_T < LF.m_T; ++jKF) {
        if (Nk > 0 && !m_marksTmp1[jKF]) {
          continue;
        }
        const Rigid3D _C = m_CsKF[jKF];
        const float imgMotion = ComputeImageMotion(z, C, _C, &first);
        if (imgMotion > imgMotionNearest) {
          continue;
        }
        imgMotionNearest = imgMotion;
        iKFNearest = jKF;
      }
      LF.m_iKFNearest = iKFNearest;
    } else {
      LF.m_iKFNearest = iKF - 1;
    }
  }
  m_Zp.DeleteKeyFrame(iKF);
  if (m_Zp.Pose::Invalid()) {
    const int iLF = m_ic2LF.front(), iKFr = m_LFs[iLF].m_iKFNearest;

    if (m_LFs[iLF].m_T.m_iFrm == m_KFs[iKFr].m_T.m_iFrm) {
      m_Zp.Initialize(m_ZpLF);
    } else {
      m_Zp.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL, iKFr, m_CsKF[iKFr],
                      BA_VARIANCE_PRIOR_GRAVITY_NEW, m_ZpLF, false, &m_CsLF[iLF].m_T,
                      BA_VARIANCE_PRIOR_POSITION_NEW, BA_VARIANCE_PRIOR_ROTATION_NEW);

    }
  }

  m_GM->LBA_DeleteKeyFrame(iFrm, iKF);
  m_GBA->PushDeleteKeyFrame(iFrm, iKF);
  //if (iKF == nKFs2) {
  //  m_usKF.Insert(0, m_usKFLast, &m_work);
  //}
#ifdef CFG_CHECK_REPROJECTION
  m_esKF.erase(m_esKF.begin() + iKF);
#endif

}

void LocalBundleAdjustor::DeleteMapPoints(const std::vector<int> &ids) {
   const int nKFs = static_cast<int>(m_KFs.size());
  std::vector<ubyte> &mcs = m_marksTmp1;
  mcs.assign(nKFs, 0);
  std::vector<ubyte> &mds = m_marksTmp2;
  mds.assign(m_ds.size(), 0);
  const int N = static_cast<int>(ids.size());

  for (int i1 = 0, i2 = 0, iKF = 0; i1 < N; i1 = i2) {
    iKF = static_cast<int>(std::upper_bound(m_iKF2d.begin() + iKF, m_iKF2d.end(), ids[i1]) -
                                            m_iKF2d.begin()) - 1;
    const int id2 = m_iKF2d[iKF + 1];

    mcs[iKF] = 1;
    for (i2 = i1 + 1; i2 < N && ids[i2] < id2; ++i2);
    for (int i = i1; i < i2; ++i) {
      const int id = ids[i];

      if (m_uds[id] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
        m_uds[id] = LBA_FLAG_TRACK_INVALID | LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      } else {
        m_uds[id] = LBA_FLAG_TRACK_INVALID;
      }
      mds[id] = 1;
    }
  }
  std::vector<int> &izsDel = m_idxsTmp1, &izs = m_idxsTmp2;
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    KeyFrame &KF = m_KFs[iKF];
    const int id = m_iKF2d[iKF], Nx = static_cast<int>(KF.m_xs.size());

    if (mcs[iKF]) {
      KF.InvalidateFeatures(mds.data() + id);
    }
    izsDel.resize(0);
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      if (!mcs[Z.m_iKF]) {
        continue;
      }
      const ubyte *mxs = mds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        if (mxs[KF.m_zs[iz].m_ix]) {
          izsDel.push_back(iz);
        }
      }
    }
    if (izsDel.empty()) {
      continue;
    }
    KF.DeleteFeatureMeasurementsPrepare(izsDel, &izs);
    KF.DeleteFeatureMeasurements(izs);
  }

  const int nLFs = static_cast<int>(m_LFs.size());
  const int STL = std::min(nLFs, LBA_MAX_SLIDING_TRACK_LENGTH);
  std::vector<std::vector<int> > &izsList = m_idxsListTmp;
  izsList.resize(STL);

  for (int ic = nLFs - 1, il = 0; ic >= 0; --ic, il = (il + STL - 1) % STL) {
    const int iLF = m_ic2LF[ic];
    LocalFrame &LF = m_LFs[iLF];
    izsDel.resize(0);
    Camera::Factor::Unitary::CC &SAczz = m_SAcusLF[iLF], &SMczz = m_SMcusLF[iLF];
    const int NZ = static_cast<int>(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      if (!mcs[Z.m_iKF]) {
        continue;
      }
      const ubyte *mxs = mds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        if (!mxs[LF.m_zs[iz].m_ix]) {
          continue;
        }
        izsDel.push_back(iz);
        Camera::Factor::Unitary::CC &Aczz = LF.m_Azs2[iz].m_Aczz;
        Aczz.MakeMinus();
        SAczz += Aczz;
        if (!(LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_NON_ZERO)) {
          continue;
        }
        Camera::Factor::Unitary::CC &Mczz = LF.m_Mzs2[iz].m_Mczz;
        Mczz.MakeMinus();
        SMczz += Mczz;
      }
    }

    std::vector<int> &izs = izsList[il];
    if (izsDel.empty()) {
      izs.resize(0);
    } else {
      LF.DeleteFeatureMeasurementsPrepare(izsDel, &izs);
      LF.DeleteFeatureMeasurements(izs);
    }
    const int NI = static_cast<int>(LF.m_Zm.m_Is.size());
    for (int iI = 0; iI < NI; ++iI) {
      const MeasurementMatchLF::Index &I = LF.m_Zm.m_Is[iI];
      //if (!mcs[I.m_iKF]) {
      //  continue;
      //}
      const int ik = I.m_ik, _il = (il + ik + 1) % STL;

      LF.m_Zm.DeleteFeatureMeasurementMatches(iI, izs, izsList[_il]);

    }

  }
  m_GBA->PushDeleteMapPoints(m_LFs[m_ic2LF.back()].m_T.m_iFrm, ids);

}

void LocalBundleAdjustor::MergeMapPoints(const std::vector<std::pair<int, int> > &ids) {
}

void LocalBundleAdjustor::UpdateCameras(const std::vector<GlobalMap::InputCamera> &Cs) {
  std::vector<int>::iterator i = m_iFrmsKF.begin();
  const int N = static_cast<int>(Cs.size()), nKFs = static_cast<int>(m_KFs.size());
  m_CsKFBkp.Resize(nKFs);
  m_marksTmp1.assign(nKFs, GM_FLAG_FRAME_DEFAULT);
  for (int j = 0; j < N; ++j) {
    const GlobalMap::InputCamera &C = Cs[j];
    i = std::lower_bound(i, m_iFrmsKF.end(), C.m_iFrm);
    if (i == m_iFrmsKF.end()) {
      break;
    } else if (*i != C.m_iFrm) {
      continue;
    }
    const int iKF = static_cast<int>(i - m_iFrmsKF.begin());
    m_CsKFBkp[iKF] = m_CsKF[iKF];
    m_CsKF[iKF] = C.m_C;
    m_marksTmp1[iKF] = GM_FLAG_FRAME_UPDATE_CAMERA;
  }
  m_GBA->PushUpdateCameras(m_LFs[m_ic2LF.back()].m_T.m_iFrm, Cs);
}

void LocalBundleAdjustor::UpdateCameras(const std::vector<ubyte> &ucs,
                                        const AlignedVector<Rigid3D> &CsKF1,
                                        const AlignedVector<Rigid3D> &CsKF2) {
  Rigid3D TI;
  int iKFNearest = -1;
  const ubyte ucmFlag = LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                        LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION |
                        LBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int ic = 0; ic < nLFs; ++ic) {
    const int iLF = m_ic2LF[ic];
    const LocalFrame &LF = m_LFs[iLF];
    const int _iKFNearest = LF.m_iKFNearest;
    if (/*_iKFNearest == -1 || */!(ucs[_iKFNearest] & GM_FLAG_FRAME_UPDATE_CAMERA)) {
      continue;
    }
    if (_iKFNearest != iKFNearest) {
      iKFNearest = _iKFNearest;
      TI = CsKF2[iKFNearest].GetInverse() * CsKF1[iKFNearest];
    }
    Camera &C = m_CsLF[iLF];
    /*if (LF.m_T.m_iFrm == m_iFrmsKF[iKFNearest]) {
      C.m_T = CsKF2[iKFNearest];
    } else*/ {
      C.m_T = C.m_T / TI;
    }
    C.m_T.GetPosition(C.m_p);
    C.m_v = TI.GetAppliedRotation(C.m_v);
    m_ucsLF[iLF] |= LBA_FLAG_FRAME_UPDATE_CAMERA;
    m_ucmsLF[iLF] |= ucmFlag;
    m_UcsLF[iLF] = LM_FLAG_FRAME_UPDATE_CAMERA_LF;
  }

  const bool ud = LBA_RESET_DEPTH_INFORMATION;
  const ubyte ucFlag = LBA_FLAG_FRAME_UPDATE_CAMERA |
                       (ud ? LBA_FLAG_FRAME_UPDATE_DEPTH : LBA_FLAG_FRAME_DEFAULT);
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(ucs[iKF] & GM_FLAG_FRAME_UPDATE_CAMERA)) {
      continue;
    }
    m_ucsKF[iKF] |= ucFlag;
    m_UcsKF[iKF] |= LM_FLAG_FRAME_UPDATE_CAMERA_KF;
    if (!ud) {
      continue;
    }
    const KeyFrame &KF = m_KFs[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      m_ucsKF[Z.m_iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
      ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        uds[KF.m_zs[iz].m_ix] |= LBA_FLAG_TRACK_UPDATE_DEPTH;
      }
    }
  }
}

void LocalBundleAdjustor::SearchMatchingKeyFrames(FRM::Frame &F) {
  std::vector<int> &iKF2Z = m_idxsTmp1;
  const int nKFs1 = static_cast<int>(m_KFs.size());
  const int nKFs2 = !F.m_Zs.empty() && F.m_Zs.back().m_iKF < nKFs1 ? nKFs1 : nKFs1 + 1;
  iKF2Z.assign(nKFs2, -1);
  const int NZ = static_cast<int>(F.m_Zs.size());
  for (int iZ = 0; iZ < NZ; ++iZ) {
    iKF2Z[F.m_Zs[iZ].m_iKF] = iZ;
  }
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = F.m_Zs[iZ];
    const GlobalMap::KeyFrame &KF = *(Z.m_iKF == nKFs1
       ? (GlobalMap::KeyFrame *) &m_IKFs2.front()
       : (GlobalMap::KeyFrame *) &m_KFs[Z.m_iKF]);
    m_marksTmp1.assign(KF.m_xs.size(), 0);
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      m_marksTmp1[F.m_zs[iz].m_ix] = 1;
    }
    const int nKFsMatch = static_cast<int>(KF.m_iKFsMatch.size());
    for (int i = 0; i < nKFsMatch; ++i) {
      const int _iKF = KF.m_iKFsMatch[i];
      if (iKF2Z[_iKF] != -1) {
        continue;
      }
      const KeyFrame &_KF = m_KFs[_iKF];
      const int _iZ = _KF.SearchFrameMeasurement(Z.m_iKF);
      if (_iZ == -1) {
        continue;
      }
      const FRM::Measurement &_Z = _KF.m_Zs[_iZ];
      const int _iz1 = _Z.m_iz1, _iz2 = _Z.m_iz2;
      int _iz;
      for (_iz = _iz1; _iz < _iz2 && !m_marksTmp1[_KF.m_zs[_iz].m_ix]; ++_iz);
      if (_iz < _iz2) {
        iKF2Z[_iKF] = -2;
      }
    }
  }
  F.m_iKFsMatch.resize(0);
  for (int iKF = 0; iKF < nKFs2; ++iKF) {
    const int iZ = iKF2Z[iKF];
    if (iZ == -1) {
      continue;
    } else if (iZ >= 0) {
      F.m_Zs[iZ].m_ik = static_cast<int>(F.m_iKFsMatch.size());
    }
    F.m_iKFsMatch.push_back(iKF);
  }
}

void LocalBundleAdjustor::PushFeatureMeasurementMatchesFirst(const FRM::Frame &F,
                                                             std::vector<int> &iKF2X, std::vector<int> &iX2z) {
  int SNx = 0;
  const int NZ = int(F.m_Zs.size());
  iKF2X.assign(m_KFs.size(), -1);
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const int iKF = F.m_Zs[iZ].m_iKF;
    iKF2X[iKF] = SNx;
    SNx += int(m_KFs[iKF].m_xs.size());
  }
  iX2z.assign(SNx, -1);
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = F.m_Zs[iZ];
    int *ix2z = iX2z.data() + iKF2X[Z.m_iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      ix2z[F.m_zs[iz].m_ix] = iz;
    }
  }
}

void LocalBundleAdjustor::PushFeatureMeasurementMatchesNext(const FRM::Frame &F1,
                                                            const FRM::Frame &F2, const std::vector<int> &iKF2X, const std::vector<int> &iX2z2,
                                                            MeasurementMatchLF &Zm) {
  ubyte firstKF = 1;
  const int NZ1 = int(F1.m_Zs.size());
  for (int iZ1 = 0; iZ1 < NZ1; ++iZ1) {
    const FRM::Measurement &Z1 = F1.m_Zs[iZ1];
    const int iX = iKF2X[Z1.m_iKF];
    if (iX == -1) {
      continue;
    }
    m_izmsTmp.resize(0);
    const int *ix2z2 = iX2z2.data() + iX;
    const int iz11 = Z1.m_iz1, iz12 = Z1.m_iz2;
    for (int iz1 = iz11; iz1 < iz12; ++iz1) {
      const int iz2 = ix2z2[F1.m_zs[iz1].m_ix];
      if (iz2 != -1) {
        m_izmsTmp.push_back(FTR::Measurement::Match(iz1, iz2));
      }
    }
    if (!m_izmsTmp.empty()) {
      Zm.PushFeatureMeasurementMatches(m_izmsTmp, Z1.m_iKF, &firstKF);
    }
  }
  if (firstKF) {
    m_izmsTmp.resize(0);
    Zm.FRM::MeasurementMatch::PushFeatureMeasurementMatches(m_izmsTmp);
  }
}

void LocalBundleAdjustor::MarkFeatureMeasurements(const LocalFrame &LF, const int iKF,
                                                  std::vector<int> &ix2z) {
  ix2z.assign(m_KFs[iKF].m_xs.size(), -1);
  const int iZ = LF.SearchFrameMeasurement(iKF);
  if (iZ == -1) {
    return;
  }
  const FRM::Measurement &Z = LF.m_Zs[iZ];
  for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
    ix2z[LF.m_zs[iz].m_ix] = iz;
  }
}

// void LocalBundleAdjustor::MarkFeatureMeasurementsUpdateDepth(const FRM::Frame &F,
//                                                              std::vector<ubyte> &ucsKF,
//                                                              std::vector<ubyte> &uds) {
//   const int NZ = static_cast<int>(F.m_Zs.size());
//   for (int iZ = 0; iZ < NZ; ++iZ) {
//     const FRM::Measurement &Z = F.m_Zs[iZ];
//     const int iKF = Z.m_iKF;
//     ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_DEPTH;
//     ubyte *uds = m_udsGT.data() + m_iKF2d[iKF];
//     for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
//       uds[F.m_zs[iz].m_ix] |= LBA_FLAG_TRACK_UPDATE_DEPTH;
//     }
//   }
// }



int LocalBundleAdjustor::CountMeasurementsFrameLF() {
  int SN = 0;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    SN += static_cast<int>(m_LFs[iLF].m_Zs.size());
  }
  return SN;
}

int LocalBundleAdjustor::CountMeasurementsFrameKF() {
  int SN = 0;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    SN += static_cast<int>(m_KFs[iKF].m_Zs.size());
  }
  return SN;
}

int LocalBundleAdjustor::CountMeasurementsFeatureLF() {
  int SN = 0;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    SN += static_cast<int>(m_LFs[iLF].m_zs.size());
  }
  return SN;
}

int LocalBundleAdjustor::CountMeasurementsFeatureKF() {
  int SN = 0;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    SN += static_cast<int>(m_KFs[iKF].m_zs.size());
  }
  return SN;
}

int LocalBundleAdjustor::CountLocalTracks() {
  int SN = 0;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (KF.m_ix2ST[ix] != KF.m_ix2ST[ix + 1]) {
        ++SN;
      }
    }
  }
  return SN;
}

int LocalBundleAdjustor::CountSlidingTracks() {
  int SN = 0;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    SN += static_cast<int>(m_KFs[iKF].m_STs.size());
  }
  return SN;
}

int LocalBundleAdjustor::CountSchurComplements() {
  return m_SAcusLF.Size() + CountSchurComplementsOffDiagonal();
}

int LocalBundleAdjustor::CountSchurComplementsOffDiagonal() {
  int SN = 0;
  const int nLFs = int(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    SN += int(m_LFs[iLF].m_iLFsMatch.size());
  }
  return SN;
}

float LocalBundleAdjustor::ComputeImageMotion(const float z1, const Rigid3D &C1, const Rigid3D &C2,
                                              ubyte *first /* = NULL */) {
  Rigid3D C12;
  C1.ApplyRotation(C2.r_20_21_22_x(), C12.r_20_21_22_x());
  if (fabs(C12.r20()) < fabs(C12.r21())) {
    C12.r_10_11_12_x().vset_all_lane(0.0f, C12.r22(), -C12.r21(), 0.0f);
    C12.r_10_11_12_x().normalize012();
    SIMD::Cross012(C12.r_10_11_12_x(), C12.r_20_21_22_x(), C12.r_00_01_02_x());
  } else {
    C12.r_00_01_02_x().vset_all_lane(C12.r22(), 0.0f, -C12.r20(), 0.0f);
    C12.r_00_01_02_x().normalize012();
    SIMD::Cross012(C12.r_20_21_22_x(), C12.r_00_01_02_x(), C12.r_10_11_12_x());
  }
  C12.SetPosition(C1.GetApplied(C2.GetPosition()));
  m_work.Resize(4 * (sizeof(Point3D) + 2 * sizeof(Point2D)) / sizeof(float));
  AlignedVector<Point3D> X1s((Point3D *) m_work.Data(), 4, false);
  AlignedVector<Point2D> x1s((Point2D *) (X1s.Data() + 4), 4, false);
  AlignedVector<Point2D> x2s(x1s.Data() + 4, 4, false);
  if (!first || *first) {
    if (first) {
      *first = 0;
    }
    const Intrinsic &K = m_K.m_K;
    x1s[0].Set(-K.fxIcx(), -K.fyIcy());
    x1s[1].Set(-K.fxIcx(), K.fyIcy());
    x1s[2].Set(K.fxIcx(), K.fyIcy());
    x1s[3].Set(K.fxIcx(), -K.fyIcy());
    for (int i = 0; i < 4; ++i) {
      X1s[i].Set(x1s[i], z1);
    }
  }
  for (int i = 0; i < 4; ++i) {
    C12.GetApplied(X1s[i]).Project(x2s[i]);
  }
  LA::AlignedVectorXf e2s((float *) x2s.Data(), 4 * 2, false);
  e2s -= (float *) x1s.Data();
  return sqrtf(e2s.SquaredLength() * m_K.m_K.fxy() * 0.25f);
}
