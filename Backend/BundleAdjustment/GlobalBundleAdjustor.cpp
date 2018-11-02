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
#include "GlobalBundleAdjustor.h"
#include "IBA_internal.h"

#ifdef CFG_DEBUG
#ifdef CFG_DEBUG_EIGEN
#define GBA_DEBUG_EIGEN
#endif
#if WIN32
//#define GBA_DEBUG_CHECK
//#define GBA_DEBUG_PRINT
#ifdef CFG_GROUND_TRUTH
//#define GBA_DEBUG_GROUND_TRUTH_MEASUREMENT
#endif
#define GBA_DEBUG_VIEW
#ifdef GBA_DEBUG_VIEW
#include "ViewerIBA.h"
#endif
#endif
#endif

#ifdef GBA_DEBUG_CHECK
//static const int g_iKF = 87;
static const int g_iKF = INT_MAX;
static const int g_ix = 44;
static const Rigid3D *g_C = NULL;
static const Camera *g_CLM = NULL;
static const GlobalBundleAdjustor::KeyFrame *g_KF = NULL;
static const Camera::Factor::Unitary::CC *g_SAcu = NULL;
static const Camera::Factor *g_SAcmLM = NULL;
static const Depth::InverseGaussian *g_d = NULL;
static const FTR::Factor::Full::Source::A *g_Ax = NULL;
static const FTR::Factor::Full::Source::M1 *g_Mx1 = NULL;
static const FTR::Factor::Full::Source::M2 *g_Mx2 = NULL;
#endif

void GlobalBundleAdjustor::Initialize(IBA::Solver *solver, const int serial, const int verbose,
                                      const int debug, const int history) {
  MT::Thread::Initialize(serial, 5, "GBA");
  m_solver = solver;
  m_LM = &solver->m_internal->m_LM;
  m_GM = &solver->m_internal->m_GM;
  m_LBA = &solver->m_internal->m_LBA;
  m_K = solver->m_internal->m_K;
  m_verbose = verbose;
  m_debug = debug;
  m_history = history;
  m_dir = solver->m_internal->m_dir;
}

void GlobalBundleAdjustor::Reset() {
  MT::Thread::Reset();
  MT_WRITE_LOCK_BEGIN(m_MT, MT_TASK_NONE, MT_TASK_GBA_Reset);
  m_ITs1.resize(0);
  m_ITs2.resize(0);
  m_IKFs1.resize(0);
  m_IKFs2.resize(0);
  m_IDKFs1.resize(0);
  m_IDKFs2.resize(0);
  m_IDMPs1.resize(0);
  m_IDMPs2.resize(0);
  m_IUCs1.resize(0);
  m_IUCs2.resize(0);
  m_IZps1.resize(0);
  m_IZps2.resize(0);
  m_IZpLM1.Invalidate();
  m_IZpLM2.Invalidate();
  MT_WRITE_LOCK_END(m_MT, MT_TASK_NONE, MT_TASK_GBA_Reset);

  for (int i = 0; i < TM_TYPES; ++i) {
    m_ts[i].Reset(TIME_AVERAGING_COUNT);
  }
  m_CsDel.resize(0);

  m_delta2 = BA_DL_RADIUS_INITIAL;

  m_Zps.resize(0);
  m_Aps.resize(0);
  //m_ZpLM.Invalidate();
  m_ZpLM.m_iKF = 0;
  m_ZpLM.Initialize(BA_WEIGHT_PRIOR_CAMERA_INITIAL,
                    BA_VARIANCE_PRIOR_VELOCITY_FIRST,
                    BA_VARIANCE_PRIOR_BIAS_ACCELERATION_FIRST,
                    BA_VARIANCE_PRIOR_BIAS_GYROSCOPE_FIRST);
  m_ApLM.MakeZero();

  m_KFs.resize(0);
  m_iFrms.resize(0);
  m_Cs.Resize(0);
  m_CsLM.Resize(0);
  m_ucs.resize(0);
  m_ucmsLM.resize(0);
#ifdef CFG_INCREMENTAL_PCG
  m_xcs.Resize(0);
  m_xmsLM.Resize(0);
#endif
  m_DsLM.Resize(0);
  m_AdsLM.Resize(0);

  m_Afps.Resize(0);
  m_AfmsLM.Resize(0);

  m_iKF2d.assign(1, 0);
  m_ds.resize(0);
  m_uds.resize(0);

  m_SAcus.Resize(0);
  m_SMcus.Resize(0);
  m_SAcmsLM.Resize(0);

  m_iKF2cb.assign(1, 0);

  //m_F = 0.0f;
}

void GlobalBundleAdjustor::PushKeyFrame(const GlobalMap::InputKeyFrame &IKF,
                                        const AlignedVector<IMU::Measurement> &us,
                                        const std::vector<Depth::InverseGaussian> &dzs
                                      ) {
  MT_WRITE_LOCK_BEGIN(m_MT, IKF.m_T.m_iFrm, MT_TASK_GBA_PushKeyFrame);
  m_ITs1.push_back(IT_KEY_FRAME);
  m_IKFs1.push_back(InputKeyFrame(IKF, us, dzs));
  MT_WRITE_LOCK_END(m_MT, IKF.m_T.m_iFrm, MT_TASK_GBA_PushKeyFrame);
}

void GlobalBundleAdjustor::PushDeleteKeyFrame(const int iFrm, const int iKF) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_PushDeleteKeyFrame);
  m_ITs1.push_back(IT_DELETE_KEY_FRAME);
  m_IDKFs1.push_back(iKF);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_PushDeleteKeyFrame);
}

void GlobalBundleAdjustor::PushDeleteMapPoints(const int iFrm, const std::vector<int> &ids) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_PushDeleteMapPoints);
  m_ITs1.push_back(IT_DELETE_MAP_POINTS);
  m_IDMPs1.push_back(ids);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_PushDeleteMapPoints);
}

void GlobalBundleAdjustor::PushUpdateCameras(const int iFrm,
                                             const std::vector<GlobalMap::InputCamera> &Cs) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_PushUpdateCameras);
  m_ITs1.push_back(IT_UPDATE_CAMERAS);
  m_IUCs1.push_back(Cs);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_PushUpdateCameras);
}

void GlobalBundleAdjustor::PushCameraPriorPose(const int iFrm, const CameraPrior::Pose &IZp) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_PushCameraPriorPose);
  m_ITs1.push_back(IT_CAMERA_PRIOR_POSE);
  m_IZps1.push_back(IZp);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_PushCameraPriorPose);
}

void GlobalBundleAdjustor::PushCameraPriorMotion(const int iFrm, const int iKF, const CameraPrior::Motion &IZp) {
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_PushCameraPriorMotion);
  std::list<InputType>::iterator IT1, IT2 = m_ITs1.begin();
  while (IT2 != m_ITs1.end()) {
    IT1 = IT2++;
    if (*IT1 == IT_CAMERA_PRIOR_MOTION) {
      m_ITs1.erase(IT1);
    }
  }
  m_ITs1.push_back(IT_CAMERA_PRIOR_MOTION);
  m_IZpLM1.Set(iKF, IZp);
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_PushCameraPriorMotion);
}

void GlobalBundleAdjustor::Run() {
  m_delta2 = BA_DL_RADIUS_INITIAL;
  m_ts[TM_TOTAL].Start();
  m_ts[TM_SYNCHRONIZE].Start();
  SynchronizeData();
  m_ts[TM_SYNCHRONIZE].Stop();
  m_ts[TM_TOTAL].Stop();

  const int iFrm = m_iFrms.back();

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
    if (GBA_EMBEDDED_POINT_ITERATION) {
      m_ts[TM_UPDATE].Start();
      EmbeddedPointIteration(m_Cs, m_ucs, m_uds, &m_ds);
      m_ts[TM_UPDATE].Stop();
    }
    MT_READ_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_BufferDataEmpty);
    m_empty = BufferDataEmpty();
    MT_READ_LOCK_END(m_MT, iFrm, MT_TASK_GBA_BufferDataEmpty);
    if (!m_empty) {
      break;
    }
  }

  m_ts[TM_TOTAL].Stop();
  for (int i = 0; i < TM_TYPES; ++i) {
    m_ts[i].Finish();
  }


  UpdateData();


  // Trigger GBA callback function if set
  if (m_callback) {
    const FRM::Tag &T = m_KFs.back().m_T;
    m_callback(T.m_iFrm, T.m_t);
  }
}

void GlobalBundleAdjustor::SetCallback(const IBA::Solver::IbaCallback& iba_callback) {
  m_callback = iba_callback;
}

//int GlobalBundleAdjustor::GetTotalPoints(int *N) {
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

float GlobalBundleAdjustor::GetTotalTime(int *N) {
  if (N) {
    *N = 0;
  }
  return 0.0f;
}

bool GlobalBundleAdjustor::SaveTimes(const std::string fileName) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  UT::PrintSaved(fileName);
  fclose(fp);
  return true;
}

bool GlobalBundleAdjustor::SaveCameras(const std::string fileName, const bool poseOnly) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  const int Nc1 = static_cast<int>(m_CsDel.size()), Nc2 = m_Cs.Size(), Nc = Nc1 + Nc2;
  AlignedVector<HistoryCamera> Cs;
  m_work.Resize(Cs.BindSize(Nc) / sizeof(float));
  Cs.Bind(m_work.Data(), Nc);
  Cs.Set(m_CsDel.data(), Nc1);
  for (int ic = 0; ic < Nc2; ++ic) {
    Cs.Push(HistoryCamera(m_KFs[ic].m_T, m_Cs[ic]));
  }
  std::sort(Cs.Data(), Cs.Data() + Nc);
  std::vector<int> &iFrm2m = m_idxsTmp1;
  if (!poseOnly) {
    iFrm2m.assign(Cs.Back().m_iFrm + 1, -1);
    for (int ic = Nc2 - m_CsLM.Size(), im = 0; ic < Nc2; ++ic, ++im) {
      iFrm2m[m_iFrms[ic]] = im;
    }
  }

  Point3D p;
  Quaternion q;
  Rotation3D R;
  LA::AlignedVector3f ba, bw;
  const Rotation3D RuT = m_K.m_Ru.GetTranspose();
  for (int ic = 0; ic < Nc2; ++ic) {
    const HistoryCamera &C = Cs[ic];
    const Rigid3D &T = C.m_C;
    T.GetPosition(p);
    p += T.GetAppliedRotationInversely(m_K.m_pu);
    Rotation3D::AB(RuT, T, R);
    R.GetQuaternion(q);
    fprintf(fp, "%f %f %f %f %f %f %f %f", C.m_t, p.x(), p.y(), p.z(),
                                           q.x(), q.y(), q.z(), q.w());
    if (!poseOnly && iFrm2m[C.m_iFrm] != -1) {
      const Camera &_C = m_CsLM[iFrm2m[C.m_iFrm]];
      RuT.Apply(_C.m_ba, ba);
      RuT.Apply(_C.m_bw, bw);
      ba += m_K.m_ba;
      bw += m_K.m_bw;
      fprintf(fp, " %f %f %f %f %f %f %f %f %f", _C.m_v.x(), _C.m_v.y(), _C.m_v.z(),
                                                 ba.x(), ba.y(), ba.z(),
                                                 bw.x(), bw.y(), bw.z());
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool GlobalBundleAdjustor::SaveCosts(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

bool GlobalBundleAdjustor::SaveResiduals(const std::string fileName, const int type) {
  FILE *fp = fopen(fileName.c_str(), "w");
  if (!fp) {
    return false;
  }
  fclose(fp);
  UT::PrintSaved(fileName);
  return true;
}

void GlobalBundleAdjustor::ComputeErrorFeature(float *ex) {
  Rigid3D Tr[2];
  FTR::Error e;
  float exi;

  *ex = 0.0f;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {

    const KeyFrame &KF = m_KFs[iKF];
    const int Nz = static_cast<int>(KF.m_zs.size());

    m_work.Resize(Nz);
    LA::AlignedVectorXf e2s(m_work.Data(), Nz, false);
    e2s.Resize(0);


    const Rigid3D C = m_Cs[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      *Tr = C / m_Cs[Z.m_iKF];

      const Depth::InverseGaussian *_ds = m_ds.data() + m_iKF2d[Z.m_iKF];
      const KeyFrame &_KF = m_KFs[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = KF.m_zs[iz];
        const int ix = z.m_ix;
        FTR::GetError(Tr, _KF.m_xs[ix], _ds[ix], z, e);

        const float e2 = e.m_e.SquaredLength();

        e2s.Push(e2);

      }
    }


    if (!e2s.Empty()) {
      const int ith = e2s.Size() >> 1;
      std::nth_element(e2s.Data(), e2s.Data() + ith, e2s.End());
      exi = sqrtf(e2s[ith] * m_K.m_K.fxy());
    } else
      exi = 0.0f;

    *ex = std::max(exi, *ex);
  }
}

void GlobalBundleAdjustor::ComputeErrorIMU(float *er, float *ep, float *ev,
                                           float *eba, float *ebw) {
  *er = *ep = *ev = *eba = *ebw = 0.0f;
  const int Nc = m_Cs.Size();
  for (int ic1 = Nc - m_CsLM.Size(), ic2 = ic1 + 1, im1 = 0, im2 = 1; ic2 < Nc; 
       ic1 = ic2++, im1 = im2++) {
    if (m_KFs[ic2].m_us.Empty()) {
      continue;
    }
    const IMU::Delta::Error e = m_DsLM[im2].GetError(m_CsLM[im1], m_CsLM[im2], m_K.m_pu,
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

void GlobalBundleAdjustor::ComputeErrorDrift(float *er, float *ep) {
  *er = 0.0f;
  *ep = 0.0f;
}

float GlobalBundleAdjustor::ComputeRMSE() {
  float Se2 = 0.0f;
  const int Nc1 = static_cast<int>(m_CsDel.size()), Nc2 = m_Cs.Size(), Nc = Nc1 + Nc2;
  return sqrtf(Se2 / Nc2);
}

void GlobalBundleAdjustor::SynchronizeData() {
  const int iFrm = m_iFrms.empty() ? MT_TASK_NONE : m_iFrms.back();
  MT_WRITE_LOCK_BEGIN(m_MT, iFrm, MT_TASK_GBA_SynchronizeData);
  m_ITs1.swap(m_ITs2);
  m_IKFs1.swap(m_IKFs2);
  m_IDKFs1.swap(m_IDKFs2);
  m_IDMPs1.swap(m_IDMPs2);
  m_IUCs1.swap(m_IUCs2);
  m_IZps1.swap(m_IZps2);
  if (m_IZpLM1.Valid()) {
    m_IZpLM2 = m_IZpLM1;
    m_IZpLM1.Invalidate();
  }
  MT_WRITE_LOCK_END(m_MT, iFrm, MT_TASK_GBA_SynchronizeData);
  m_Ucs.assign(m_KFs.size(), GM_FLAG_FRAME_DEFAULT);

  while (!m_ITs2.empty()) {
    const InputType IT = m_ITs2.front();
    m_ITs2.pop_front();
    if (IT == IT_KEY_FRAME) {
      PushKeyFrame(m_IKFs2.front());
      m_IKFs2.pop_front();
    } else if (IT == IT_DELETE_KEY_FRAME) {
      DeleteKeyFrame(m_IDKFs2.front());
      m_IDKFs2.pop_front();
    } else if (IT == IT_DELETE_MAP_POINTS) {
      DeleteMapPoints(m_IDMPs2.front());
      m_IDMPs2.pop_front();
    } else if (IT == IT_UPDATE_CAMERAS) {
      UpdateCameras(m_IUCs2.front());
      m_IUCs2.pop_front();
    } else if (IT == IT_CAMERA_PRIOR_POSE) {
      PushCameraPriorPose(m_IZps2.front());
      m_IZps2.pop_front();
    } else if (IT == IT_CAMERA_PRIOR_MOTION) {
      if (m_IZpLM2.Valid()) {
        SetCameraPriorMotion(m_IZpLM2);
        m_IZpLM2.Invalidate();
      }
    }
  }
}

void GlobalBundleAdjustor::UpdateData() {
  m_GM->GBA_Update(m_iFrms, m_Cs, m_Ucs);
}

bool GlobalBundleAdjustor::BufferDataEmpty() {
  return m_IKFs1.empty() && m_IDKFs1.empty() && m_IUCs1.empty();
}

void GlobalBundleAdjustor::PushKeyFrame(const InputKeyFrame &IKF) {
  //Timer timer;
  //timer.Start();
  const int nKFs1 = static_cast<int>(m_KFs.size()), nKFs2 = nKFs1 + 1;
  const int Nm1 = m_CsLM.Size(), Nm2 = Nm1 + 1;
  m_KFs.resize(nKFs2);
  m_iFrms.resize(nKFs2);
  m_Cs.Resize(nKFs2, true);
  m_CsLM.Resize(Nm2, true);

  m_ucs.resize(nKFs2, GBA_FLAG_FRAME_UPDATE_CAMERA);
  m_Ucs.resize(nKFs2, GM_FLAG_FRAME_DEFAULT);
#ifdef CFG_INCREMENTAL_PCG
  m_xcs.Resize(nKFs2, true);
  m_xmsLM.Resize(Nm2, true);
#endif
  const ubyte ucmFlag1 = GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                         GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  const ubyte ucmFlag2 = GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                         GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                         GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
  m_ucmsLM.resize(Nm2, ucmFlag1 | ucmFlag2);
  m_DsLM.Resize(Nm2, true);

  const ubyte udFlag1 = GBA_FLAG_TRACK_UPDATE_DEPTH | GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
  const ubyte udFlag2 = GBA_FLAG_TRACK_UPDATE_DEPTH
#ifdef GBA_FLAG_TRACK_MEASURE
                      | GBA_FLAG_TRACK_MEASURE
#endif
                      ;
  const int iKF = nKFs1, im = Nm1;
  KeyFrame &KF = m_KFs[iKF];
  KF.Initialize(IKF);
  const int iFrm = KF.m_T.m_iFrm;
  m_iFrms[iKF] = iFrm;
  std::vector<int> &iKF2X = m_idxsTmp1, &iX2z = m_idxsTmp2;
  PushFeatureMeasurementMatchesFirst(KF, iKF2X, iX2z);
  const int Nk = static_cast<int>(KF.m_iKFsMatch.size());

  for (int ik = 0; ik < Nk; ++ik) {
    KeyFrame &_KF = m_KFs[KF.m_iKFsMatch[ik]];
    PushFeatureMeasurementMatchesNext(_KF, KF, iKF2X, iX2z, KF.m_Zm);

    _KF.m_iKFsMatch.push_back(iKF);
  }
  if (iKF > 0) {
    const int _iKF = iKF - 1;
    if (KF.m_us.Empty()) {
      m_ucmsLM[im] |= GBA_FLAG_CAMERA_MOTION_INVALID;
    } else {
      if ((KF.m_T.m_t - m_KFs[_iKF].m_T.m_t) > GBA_MAX_IMU_TIME_INTERVAL) {
        KF.m_us.Clear();
        m_ucmsLM[im] |= GBA_FLAG_CAMERA_MOTION_INVALID;
      } else if (KF.m_iKFsMatch.empty() || KF.m_iKFsMatch.back() != _iKF) {
        KF.InsertMatchKeyFrame(_iKF);
        m_KFs[_iKF].m_iKFsMatch.push_back(iKF);
      }
      if (!KF.m_us.Empty()) {
        const int _im = im - 1;
        if (_im >= 0 && (m_ucmsLM[_im] & GBA_FLAG_CAMERA_MOTION_INVALID)) {
          m_ucmsLM[_im] &= ~GBA_FLAG_CAMERA_MOTION_INVALID;
          m_ucmsLM[_im] |= ucmFlag2;
        }
      }
    }
  }
  const bool v1 = IKF.m_C.m_T.Valid(), v2 = IKF.m_C.m_v.Valid();
  if (im > 0 && !KF.m_us.Empty() && (GBA_PROPAGATE_CAMERA || !v1 || !v2)) {
    IMU::Delta D;
    Camera C;
    const KeyFrame &_KF = m_KFs[iKF - 1];
    const Camera &_C = m_CsLM[im - 1];
    IMU::PreIntegrate(KF.m_us, _KF.m_T.m_t, KF.m_T.m_t, _C, &D, &m_work, false,
                      _KF.m_us.Empty() ? NULL : &_KF.m_us.Back(), NULL, BA_ANGLE_EPSILON);
    IMU::Propagate(m_K.m_pu, D, _C, C, BA_ANGLE_EPSILON);
    if (GBA_PROPAGATE_CAMERA) {
      m_CsLM[im] = C;
      m_Cs[iKF] = m_CsLM[im].m_T;
    } else {
      if (!v1) {
        m_Cs[iKF] = C.m_T;
        m_CsLM[im].m_T = C.m_T;
        m_CsLM[im].m_p = C.m_p;
      }
      if (!v2) {
        m_CsLM[im].m_v = C.m_v;
        m_CsLM[im].m_ba = C.m_ba;
        m_CsLM[im].m_bw = C.m_bw;
      }
    }
  } else {
    m_Cs[iKF] = IKF.m_C.m_T;
    m_CsLM[im] = IKF.m_C;
  }

#ifdef CFG_INCREMENTAL_PCG
  m_xcs[iKF].MakeZero();
  m_xmsLM[im].MakeZero();
#endif
  IMU::Delta &D = m_DsLM[im];
  if (KF.m_us.Empty() || im == 0) {
    D.Invalidate();
  } else {
    const int _iKF = iKF - 1;
    const int _im = im - 1;
    const KeyFrame &_KF = m_KFs[_iKF];
    const float _t = _KF.m_T.m_t;
    IMU::PreIntegrate(KF.m_us, _t, KF.m_T.m_t, m_CsLM[_im], &D, &m_work, true,
                      _KF.m_us.Empty() ? NULL : &_KF.m_us.Back(), NULL, BA_ANGLE_EPSILON);
    if (_im > 0 && !_KF.m_us.Empty()) {
      IMU::Delta &_D = m_DsLM[_im];
      IMU::PreIntegrate(_KF.m_us, m_KFs[_iKF - 1].m_T.m_t, _t, m_CsLM[_im - 1], &_D, &m_work, true,
                        _D.m_u1.Valid() ? &_D.m_u1 : NULL, &KF.m_us.Front(), BA_ANGLE_EPSILON);
      m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_CAMERA;
      m_ucmsLM[_im] |= ucmFlag1 | ucmFlag2;
    }
  }
  if (iFrm == 0) {
    const LA::Vector3f s2r = LA::Vector3f::Get(BA_VARIANCE_FIX_ORIGIN_ROTATION_X,
                                               BA_VARIANCE_FIX_ORIGIN_ROTATION_Y,
                                               BA_VARIANCE_FIX_ORIGIN_ROTATION_Z);
    const float s2p = BA_VARIANCE_FIX_ORIGIN_POSITION;
    m_Zo.Set(BA_WEIGHT_FIX_ORIGIN, s2r, s2p, m_Cs[iKF]);
    m_Ao.MakeZero();
  }

  m_iKF2d.push_back(m_iKF2d.back());
  m_iKF2cb.push_back(m_iKF2cb.back() + static_cast<int>(KF.m_ikp2KF.size()));
  const int NZ = static_cast<int>(IKF.m_Zs.size());
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const FRM::Measurement &Z = IKF.m_Zs[iZ];
    const int _iKF = Z.m_iKF, id = m_iKF2d[_iKF];
    Depth::InverseGaussian *ds = m_ds.data() + id;
    ubyte *uds = m_uds.data() + id/*, *Uds = m_Uds.data() + id*/;

    bool ud = false;
    const KeyFrame &_KF = m_KFs[_iKF];
    for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
      const int ix = IKF.m_zs[iz].m_ix;
      if (IKF.m_dzs[iz].Valid() && fabs(IKF.m_dzs[iz].u() - ds[ix].u()) >= BA_UPDATE_DEPTH) {
        ds[ix] = IKF.m_dzs[iz];
        m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_DEPTH;
        uds[ix] |= udFlag2;
        //m_Ucs[_iKF] |= GM_FLAG_FRAME_UPDATE_DEPTH;
        //Uds[ix] |= GM_FLAG_TRACK_UPDATE_DEPTH;
        ud = true;
      }
#ifdef GBA_FLAG_TRACK_MEASURE
      else {
        uds[ix] |= GBA_FLAG_TRACK_MEASURE;
      }
#endif

    }

  }
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
    m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_DEPTH;

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
#ifdef GBA_FLAG_TRACK_MEASURE
      if (Nz > 0) {
        uds[i] |= GBA_FLAG_TRACK_MEASURE;
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
    const int Nx1 = static_cast<int>(m_KFs[_iKF].m_xs.size());
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
    std::vector<int> &izs = m_idxsTmp2, &ix2z = m_idxsTmp3;
    izs.resize(Nk);
    KeyFrame &_KF = m_KFs[_iKF];
    for (int ik2 = 0; ik2 < Nk; ++ik2) {
      const int iKF2 = ik2KF[ik2];
      KeyFrame &KF2 = m_KFs[iKF2];
      int jk2, &iz2 = izs[ik2];
      const std::vector<FTR::Measurement> &zs2 = m_zsListTmp[ik2];
      const int Nk21 = static_cast<int>(KF2.m_iKFsMatch.size());
      KF2.PushFeatureMeasurements(_iKF, zs2, &jk2, &iz2, &m_work);
      const int Nk22 = static_cast<int>(KF2.m_iKFsMatch.size());
      if (Nk22 > Nk21) {
        const std::vector<int>::iterator _ik = std::lower_bound(_KF.m_iKFsMatch.begin() +
                                                                _KF.m_Zm.m_SMczms.Size(),
                                                                _KF.m_iKFsMatch.end(), iKF2);

        _KF.m_iKFsMatch.insert(_ik, iKF2);
      }
      const int Nz2 = static_cast<int>(zs2.size());
      KF2.m_Zm.InsertFeatureMeasurement2(jk2, iz2, Nz2);
      for (int jk2 = KF2.m_Zm.m_SMczms.Size(); jk2 < Nk22; ++jk2) {
        KeyFrame &KF3 = m_KFs[KF2.m_iKFsMatch[jk2]];
        const int jk3 = KF3.SearchMatchKeyFrame(iKF2);

        KF3.m_Zm.InsertFeatureMeasurement1(jk3, iz2, Nz2);
      }
      ix2z.assign(Nx, -1);
      int *_ix2z = ix2z.data() - Nx1;
      for (int i = 0; i < Nz2; ++i) {
        const int ix = zs2[i].m_ix;

        _ix2z[ix] = iz2 + i;
      }
      for (int ik1 = 0; ik1 < ik2; ++ik1) {
        const int iKF1 = ik2KF[ik1];
        m_izmsTmp.resize(0);
        const int iz1 = izs[ik1];
        const std::vector<FTR::Measurement> &zs1 = m_zsListTmp[ik1];
        const int Nz1 = static_cast<int>(zs1.size());
        for (int i = 0; i < Nz1; ++i) {
          const int ix = zs1[i].m_ix, _iz2 = _ix2z[ix];
          if (_iz2 == -1) {
            continue;
          }

          m_izmsTmp.push_back(FTR::Measurement::Match(iz1 + i, _iz2));
        }
        if (m_izmsTmp.empty()) {
          continue;
        }
        const std::vector<int>::iterator _jk2 = std::lower_bound(KF2.m_iKFsMatch.begin() + jk2,
                                                                 KF2.m_iKFsMatch.begin() +
                                                                 KF2.m_Zm.m_SMczms.Size(), iKF1);
        jk2 = static_cast<int>(_jk2 - KF2.m_iKFsMatch.begin());
        if (_jk2 == KF2.m_iKFsMatch.end() || *_jk2 != iKF1) {
          KF2.InsertMatchKeyFrame(iKF1, &_jk2);
          KeyFrame &KF1 = m_KFs[iKF1];
          const std::vector<int>::iterator jk1 = std::lower_bound(KF1.m_iKFsMatch.begin() +
                                                                  KF1.m_Zm.m_SMczms.Size(),
                                                                  KF1.m_iKFsMatch.end(), iKF2);

          KF1.m_iKFsMatch.insert(jk1, iKF2);
        }
        KF2.m_Zm.InsertFeatureMeasurementMatches(jk2, m_izmsTmp, &m_work);

      }
      const int Nk23 = static_cast<int>(KF2.m_iKFsMatch.size());
      if (Nk23 == Nk21) {
        continue;
      }
      const int Ncb = Nk23 - Nk21;
      for (int jKF = iKF2; jKF <= iKF; ++jKF) {
        m_iKF2cb[jKF + 1] += Ncb;
      }
    }
    m_KFs[_iKF].PushFeatures(m_xsTmp);
  }
  m_AdsLM.InsertZero(Nm1);
  m_Afps.InsertZero(nKFs1);
  m_AfmsLM.InsertZero(Nm1);
  m_SAcus.InsertZero(nKFs1);
  m_SMcus.InsertZero(nKFs1);
  m_SAcmsLM.InsertZero(Nm1);
}

void GlobalBundleAdjustor::DeleteKeyFrame(const int iKF) {
  //Timer timer;
  //timer.Start();
  const int iFrm = m_iFrms[iKF];
  const int nKFs = static_cast<int>(m_KFs.size());

  m_KFs[iKF + 1].m_us.Insert(0, m_KFs[iKF].m_us, &m_work);
  if (iKF > 0 && !m_KFs[iKF + 1].m_us.Empty() &&
      m_KFs[iKF + 1].SearchMatchKeyFrame(iKF - 1) == -1) {
    m_KFs[iKF + 1].InsertMatchKeyFrame(iKF - 1);
    KeyFrame &KF = m_KFs[iKF - 1];
    const std::vector<int>::iterator ik = std::lower_bound(KF.m_iKFsMatch.begin() +
                                                           KF.m_Zm.m_SMczms.Size(),
                                                           KF.m_iKFsMatch.end(), iKF + 1);

    KF.m_iKFsMatch.insert(ik, iKF + 1);
  }
  const int Nd = static_cast<int>(m_KFs[iKF].m_xs.size());
  const int id1 = m_iKF2d[iKF], id2 = id1 + Nd;
  const ubyte *uds = m_uds.data() + id1;
  m_iKF2Z.resize(nKFs);
  for (int jKF2 = iKF + 1; jKF2 < nKFs; ++jKF2) {
    KeyFrame &KF2 = m_KFs[jKF2];
    const std::vector<FRM::Measurement>::iterator iZ2 = std::lower_bound(KF2.m_Zs.begin(),
                                                                         KF2.m_Zs.end(), iKF);
    m_iKF2Z[jKF2] = iZ2;
    const bool z2 = iZ2 != KF2.m_Zs.end() && iZ2->m_iKF == iKF;
    //const int Nk = static_cast<int>(KF2.m_iKFsMatch.size());
    const int Nk = KF2.m_Zm.m_SMczms.Size();
    for (int ik = 0; ik < Nk; ++ik) {
      const int jKF1 = KF2.m_iKFsMatch[ik];
      const std::vector<FRM::Measurement>::iterator iZ1 = m_iKF2Z[jKF1];
      const KeyFrame &KF1 = m_KFs[jKF1];
      const bool z1 = jKF1 > iKF && iZ1 != KF1.m_Zs.end() && iZ1->m_iKF == iKF;
      if (!z1 && !z2) {
        continue;
      }
      std::vector<FTR::Measurement::Match> &izms = KF2.m_Zm.m_izms;
      const int i1 = KF2.m_Zm.m_ik2zm[ik], i2 = KF2.m_Zm.m_ik2zm[ik + 1];
      if (z1) {
        const int Nz1 = iZ1->CountFeatureMeasurements();
        for (int i = i2 - 1; i >= i1 && izms[i].m_iz1 >= iZ1->m_iz2; --i) {
          izms[i].m_iz1 -= Nz1;
        }
      }
      if (z2) {
        const int Nz2 = iZ2->CountFeatureMeasurements();
        int i;
        for (i = i2 - 1; i >= i1 && izms[i].m_iz2 >= iZ2->m_iz2; --i) {
          izms[i].m_iz2 -= Nz2;
        }
        if (!z1) {
          continue;
        }
        const int j2 = i + 1;
        const int j1 = static_cast<int>(std::lower_bound(izms.begin() + i1,
                                                         izms.begin() + j2, iZ2->m_iz1) -
                                                         izms.begin());
        const int Nzm = j2 - j1;
        m_marksTmp1.assign(Nzm, 0);
        ubyte *ms = m_marksTmp1.data() - j1;
        for (int j = j1; j < j2; ++j) {
          const int ix = KF2.m_zs[izms[j].m_iz2].m_ix;
          ms[j] = (uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) ? 0 : 1;
        }
        KF2.m_Zm.DeleteFeatureMeasurementMatches(ik, j1, j2, ms);
      }
    }
  }
  {
    KeyFrame &KF = m_KFs[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iKF = Z.m_iKF;
      Camera::Factor::Unitary::CC &SAcxx = m_SAcus[_iKF];
      ubyte *_uds = m_uds.data() + m_iKF2d[_iKF];
      KeyFrame &_KF = m_KFs[_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = KF.m_zs[iz].m_ix;
        FTR::Factor::Full::A2 &A = KF.m_Azs2[iz];
        A.m_Acxx.MakeMinus();
        SAcxx += A.m_Acxx;
        A.m_adx.MakeMinus();
        _KF.m_Axs[ix].m_Sadx += A.m_adx;
        _uds[ix] |= GBA_FLAG_TRACK_UPDATE_INFORMATION;
      }
      if (Z.m_iz1 < Z.m_iz2) {
        m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
      }
    }
  }
  int icb = m_iKF2cb[iKF];
  for (int jKF = 0; jKF < iKF; ++jKF) {
    m_KFs[jKF].DeleteMatchKeyFrame(iKF);
  }
  for (int jKF = iKF + 1; jKF < nKFs; ++jKF) {
    const std::vector<FRM::Measurement>::iterator iZ = m_iKF2Z[jKF];
    KeyFrame &KF = m_KFs[jKF];
    const int _iZ = iZ != KF.m_Zs.end() && iZ->m_iKF == iKF ?
                      static_cast<int>(iZ - KF.m_Zs.begin()) : -1;
    if (_iZ != -1) {
      Camera::Factor::Unitary::CC &SAczz = m_SAcus[jKF], &SMczz = m_SMcus[jKF];
      for (int iz = iZ->m_iz1; iz < iZ->m_iz2; ++iz) {
        Camera::Factor::Unitary::CC &Aczz = KF.m_Azs2[iz].m_Aczz;
        Aczz.MakeMinus();
        SAczz += Aczz;
        if (uds[KF.m_zs[iz].m_ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
          continue;
        }
        Camera::Factor::Unitary::CC &Mczz = KF.m_Mzs2[iz].m_Mczz;
        Mczz.MakeMinus();
        SMczz += Mczz;
      }
    }
    KF.DeleteKeyFrame(iKF, &iZ);
    m_iKF2d[jKF] -= Nd;
    m_iKF2cb[jKF] = icb;
    icb += static_cast<int>(KF.m_ikp2KF.size());
  }
  m_CsDel.push_back(HistoryCamera(m_KFs[iKF].m_T, m_Cs[iKF]));
  m_KFs.erase(m_KFs.begin() + iKF);
  m_iFrms.erase(m_iFrms.begin() + iKF);
  m_iKF2d.back() -= Nd;
  m_iKF2d.erase(m_iKF2d.begin() + iKF);
  m_iKF2cb.back() = icb;
  m_iKF2cb.erase(m_iKF2cb.begin() + iKF);
  m_Cs.Erase(iKF);
  const int Nm = m_CsLM.Size(), im = iKF - nKFs + Nm;
  if (im >= 0) {
    m_CsLM.Erase(im);
  }

#ifdef CFG_INCREMENTAL_PCG
  m_xcs.Erase(iKF);
  if (im >= 0) {
    m_xmsLM.Erase(im);
  }
#endif
  m_ucs.erase(m_ucs.begin() + iKF);
  if (im >= 0) {
    m_ucmsLM.erase(m_ucmsLM.begin() + im);
    const int im1 = im - 1, im2 = im + 1;
    if (im1 >= 0 && im2 < Nm) {
      KeyFrame &KF1 = m_KFs[iKF - 1], &KF2 = m_KFs[iKF];
      IMU::Delta &D = m_DsLM[im2];
      if ((KF2.m_T.m_t - KF1.m_T.m_t) > GBA_MAX_IMU_TIME_INTERVAL) {
        KF2.m_us.Clear();
        D.Invalidate();
        m_ucmsLM[im] |= GBA_FLAG_CAMERA_MOTION_INVALID;
      } else {
        IMU::PreIntegrate(KF2.m_us, KF1.m_T.m_t, KF2.m_T.m_t, m_CsLM[im1], &D, &m_work, false,
                          KF1.m_us.Empty() ? NULL : &KF1.m_us.Back(), NULL, BA_ANGLE_EPSILON);
      }
    }
    m_DsLM.Erase(im);
  }
  m_Ucs.erase(m_Ucs.begin() + iKF);
  m_ds.erase(m_ds.begin() + id1, m_ds.begin() + id2);
  m_uds.erase(m_uds.begin() + id1, m_uds.begin() + id2);

  m_Afps.Erase(iKF);
  m_SAcus.Erase(iKF);
  m_SMcus.Erase(iKF);
  if (im >= 0) {
    const int im1 = im - 1, im2 = im + 1;
    //Camera::Factor &SAcm = m_SAcmsLM[im];
    if (im1 >= 0) {
      IMU::Delta::Factor &A = m_AdsLM[im];
      A.m_A11.MakeMinus();
      m_SAcus[iKF - 1] += A.m_A11.m_Acc;
      Camera::Factor::Unitary &SAcm1 = m_SAcmsLM[im1].m_Au;
      SAcm1.m_Acm += A.m_A11.m_Acm;
      SAcm1.m_Amm += A.m_A11.m_Amm;
    }
    if (im2 < Nm) {
      IMU::Delta::Factor &A = m_AdsLM[im2];
      A.m_A22.MakeMinus();
      m_SAcus[iKF] += A.m_A22.m_Acc;
      Camera::Factor::Unitary &SAcm2 = m_SAcmsLM[im2].m_Au;
      SAcm2.m_Acm += A.m_A22.m_Acm;
      SAcm2.m_Amm += A.m_A22.m_Amm;
    }
    if (im1 >= 0 && im2 < Nm && m_DsLM[im].Valid()) {
      IMU::Delta::Factor::Auxiliary::Global U;
      Camera::Factor::Unitary &SAcm1 = m_SAcmsLM[im1].m_Au;
      Camera::Factor &SAcm2 = m_SAcmsLM[im2];
      IMU::Delta::Factor &A = m_AdsLM[im2];
      m_DsLM[im].GetFactor(BA_WEIGHT_IMU, m_CsLM[im1], m_CsLM[im], m_K.m_pu,
                           &A, &SAcm2.m_Ab, &U, BA_ANGLE_EPSILON);
      m_SAcus[iKF - 1] += A.m_A11.m_Acc;
      SAcm1.m_Acm += A.m_A11.m_Acm;
      SAcm1.m_Amm += A.m_A11.m_Amm;
      m_SAcus[iKF] += A.m_A22.m_Acc;
      SAcm2.m_Au.m_Acm += A.m_A22.m_Acm;
      SAcm2.m_Au.m_Amm += A.m_A22.m_Amm;
    } else if (im1 == -1) {
      m_SAcmsLM[im2].m_Ab.MakeZero();
    }
    m_AdsLM.Erase(im);
    m_AfmsLM.Erase(im);
    m_SAcmsLM.Erase(im);
  }
  int iZp, jZp;
  Camera::Factor::Unitary::CC Au, dAu;
  Camera::Factor::Binary::CC Ab, dAb;
  CameraPrior::Pose::Factor::Auxiliary Up;
  const int NZp1 = static_cast<int>(m_Zps.size());
  for (iZp = 0; iZp < NZp1; ++iZp) {
    CameraPrior::Pose &Zp = m_Zps[iZp];
    const std::vector<int>::iterator i = std::lower_bound(Zp.m_iKFs.begin(),
                                                          Zp.m_iKFs.end(), iKF);
    const int _i = static_cast<int>(i - Zp.m_iKFs.begin()) + 1;
    const bool del = Zp.m_iKFr == iKF || (i != Zp.m_iKFs.end() && *i == iKF);
    Zp.DeleteKeyFrame(iKF, &i);
    const int N = static_cast<int>(Zp.m_iKFs.size());
    if (!del && N != 0) {
      continue;
    }
    CameraPrior::Pose::Factor &Ap = m_Aps[iZp];
    Ap.Swap(m_ApTmp, m_bpTmp);
    const bool nonZero = Zp.Valid() && N != 0;
    if (nonZero) {
      m_ApTmp.Erase(_i);
      m_bpTmp.Erase(_i);
      m_work.Resize(Up.BindSize(N) / sizeof(float));
      Up.Bind(m_work.Data(), N);
      Zp.GetFactor(BA_WEIGHT_PRIOR_CAMERA_POSE, m_Cs, &Ap, &Up, BA_ANGLE_EPSILON);
    }
    for (int j = -1, _j = 0; j < N; ++j, ++_j) {
      const int jKF = j == -1 ? Zp.m_iKFr : Zp.m_iKFs[j];
      if (jKF != -1) {
        dAu.Set(m_ApTmp[_j][_j], m_bpTmp[_j]);
        if (nonZero) {
          Au.Set(Ap.m_A[_j][_j], Ap.m_b[_j]);
          Camera::Factor::Unitary::CC::AmB(Au, dAu, dAu);
        } else {
          dAu.MakeMinus();
        }
        m_SAcus[jKF] += dAu;
      }
      if (j == -1) {
        continue;
      }
      if (Zp.Valid()) {
        const int iKF1 = std::min(Zp.m_iKFr, jKF), iKF2 = std::max(Zp.m_iKFr, jKF);
        KeyFrame &KF2 = m_KFs[iKF2];
        m_ApTmp[0][_j].Get(dAb, Zp.m_iKFr > jKF);
        if (nonZero) {
          Ap.m_A[0][_j].Get(Ab, Zp.m_iKFr > jKF);
          Camera::Factor::Binary::CC::AmB(Ab, dAb, dAb);
        } else {
          dAb.MakeMinus();
        }
        const int ip = KF2.SearchPriorKeyFrame(iKF1);
        KF2.m_SAps[ip] += dAb;
      }
    }
    for (int i1 = 0, _i1 = 1; i1 < N; ++i1, ++_i1) {
      const int iKF1 = Zp.m_iKFs[i1];
      for (int i2 = i1 + 1, _i2 = i2 + 1; i2 < N; ++i2, ++_i2) {
        const int iKF2 = Zp.m_iKFs[i2];
        KeyFrame &KF = m_KFs[iKF2];
        const int ip = KF.SearchPriorKeyFrame(iKF1);
        if (nonZero) {
          Camera::Factor::Binary::CC::AmB(Ap.m_A[_i1][_i2], m_ApTmp[_i1][_i2], dAb);
        } else {
          m_ApTmp[_i1][_i2].GetMinus(dAb);
        }
        KF.m_SAps[ip] += dAb;
      }
    }
  }
  for (iZp = jZp = 0; iZp < NZp1; ++iZp) {
    const CameraPrior::Pose &Zp = m_Zps[iZp];
    if (Zp.Invalid() || Zp.m_iKFs.empty()) {
      continue;
    } else if (jZp < iZp) {
      m_Zps[jZp] = Zp;
      m_Aps[jZp] = m_Aps[iZp];
    }
    ++jZp;
  }
  const int NZp2 = jZp;
  m_Zps.resize(NZp2);
  m_Aps.resize(NZp2);
  if (m_ZpLM.DeleteKeyFrame(iKF)) {
    m_ApLM.MakeZero();
  }

}

void GlobalBundleAdjustor::DeleteMapPoints(const std::vector<int> &ids) {
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

      if (m_uds[id] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
        m_uds[id] = GBA_FLAG_TRACK_INVALID | GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
        mds[id] = 2;
      } else {
        m_uds[id] = GBA_FLAG_TRACK_INVALID;
        mds[id] = 1;
      }
    }
  }
  std::vector<int> &izsDel = m_idxsTmp1;
  std::vector<ubyte> &mzs = m_marksTmp3;
  std::vector<std::vector<int> > &izsList = m_idxsListTmp;
  izsList.resize(nKFs);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    KeyFrame &KF = m_KFs[iKF];
    const int id = m_iKF2d[iKF], Nx = static_cast<int>(KF.m_xs.size());

    if (mcs[iKF]) {
      Camera::Factor::Unitary::CC &SMcxx = m_SMcus[iKF];
      const ubyte *mxs = mds.data() + id;
      for (int ix = 0; ix < Nx; ++ix) {
        if (mxs[ix] != 1) {
          continue;
        }
        Camera::Factor::Unitary::CC &Mcxx = KF.m_Mxs2[ix].m_Mcxx;
        Mcxx.MakeMinus();
        SMcxx += Mcxx;
      }
      KF.InvalidateFeatures(mxs);
    }
    izsDel.resize(0);
    const int Nz = static_cast<int>(KF.m_zs.size());
    mzs.resize(Nz, 0);
    Camera::Factor::Unitary::CC &SAczz = m_SAcus[iKF], &SMczz = m_SMcus[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      if (!mcs[Z.m_iKF]) {
        continue;
      }
      Camera::Factor::Unitary::CC &SAcxx = m_SAcus[Z.m_iKF];
      Camera::Factor::Binary::CC &SAcxz = KF.m_SAcxzs[iZ], &SMcxz = KF.m_Zm.m_SMczms[Z.m_ik];
      const ubyte *mxs = mds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = KF.m_zs[iz].m_ix;
        if (!mxs[ix]) {
          continue;
        }
        izsDel.push_back(iz);
        mzs[iz] = mxs[ix] == 2;
        FTR::Factor::Full::A2 &Az = KF.m_Azs2[iz];
        Az.m_Acxx.MakeMinus();
        SAcxx += Az.m_Acxx;
        Az.m_Acxz.MakeMinus();
        SAcxz += Az.m_Acxz;
        Az.m_Aczz.MakeMinus();
        SAczz += Az.m_Aczz;
        if (mzs[iz]) {
          continue;
        }
        FTR::Factor::Full::M2 &Mz = KF.m_Mzs2[iz];
        Mz.m_Mcxz.MakeMinus();
        SMcxz += Mz.m_Mcxz;
        Mz.m_Mczz.MakeMinus();
        SMczz += Mz.m_Mczz;
      }
    }
    std::vector<int> &izs = izsList[iKF];
    if (izsDel.empty()) {
      izs.resize(0);
    } else {
      KF.DeleteFeatureMeasurementsPrepare(izsDel, &izs);
      KF.DeleteFeatureMeasurements(izs);
      for (int iz = 0; iz < Nz; ++iz) {
        if (izs[iz] == -1 && mzs[iz]) {
          izs[iz] = -2;
        }
      }
    }
    //const int Nk = static_cast<int>(KF.m_iKFsMatch.size());
    const int Nk = KF.m_Zm.m_SMczms.Size();
    for (int ik = 0; ik < Nk; ++ik) {
      KF.m_Zm.DeleteFeatureMeasurementMatches(ik, KF.m_Zm.m_ik2zm[ik], KF.m_Zm.m_ik2zm[ik + 1],
                                              izsList[KF.m_iKFsMatch[ik]], izs);
    }
  }
}

void GlobalBundleAdjustor::UpdateCameras(const std::vector<GlobalMap::InputCamera> &Cs) {
  std::vector<int>::iterator i = m_iFrms.begin();
  const int N = static_cast<int>(Cs.size()), nKFs = static_cast<int>(m_KFs.size());
  for (int j = 0; j < N; ++j) {
    const GlobalMap::InputCamera &C = Cs[j];
    i = std::lower_bound(i, m_iFrms.end(), C.m_iFrm);
    if (i == m_iFrms.end()) {
      break;
    } else if (*i != C.m_iFrm) {
      continue;
    }
    const int iKF = static_cast<int>(i - m_iFrms.begin());
    m_Cs[iKF] = C.m_C;
    m_ucs[iKF] |= GBA_FLAG_FRAME_UPDATE_CAMERA;
    m_Ucs[iKF] |= GM_FLAG_FRAME_UPDATE_CAMERA;
    const int im = iKF - (nKFs - m_CsLM.Size());
    if (im < 0) {
      continue;
    }
    Camera &_C = m_CsLM[im];
    _C.m_T = C.m_C;
    _C.m_T.GetPosition(_C.m_p);
    m_ucmsLM[im] |= GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                    GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  }
}

void GlobalBundleAdjustor::PushCameraPriorPose(const CameraPrior::Pose &Zp) {
  const int iZp = static_cast<int>(m_Zps.size());
  m_Zps.push_back(Zp);
  const int N = static_cast<int>(Zp.m_iKFs.size());
  m_Aps.resize(iZp + 1);
  //m_Aps.push_back(CameraPrior::Pose::Factor());
  m_Aps[iZp].MakeZero(N);
  m_ucs[Zp.m_iKFr] |= GBA_FLAG_FRAME_UPDATE_CAMERA;
  const int nKFs = static_cast<int>(m_KFs.size()), imr = Zp.m_iKFr - (nKFs - m_CsLM.Size());
  if (imr >= 0) {
    m_ucmsLM[imr] |= GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION;
  }
  const ubyte ucmFlag = GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  for (int i = 0; i < N; ++i) {
    const int _iKF = Zp.m_iKFs[i];
    m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_CAMERA;
    const int im = _iKF - (nKFs - m_CsLM.Size());
    if (im >= 0) {
      m_ucmsLM[im] |= ucmFlag;
    }
  }
  for (int i1 = -1; i1 < N; ++i1) {
    const int iKF1 = i1 == -1 ? Zp.m_iKFr : Zp.m_iKFs[i1];
    for (int i2 = i1 + 1; i2 < N; ++i2) {
      const int iKF2 = Zp.m_iKFs[i2];
      const int _iKF1 = std::min(iKF1, iKF2), _iKF2 = std::max(iKF1, iKF2);
      KeyFrame &KF = m_KFs[_iKF2];
      const int Nkp = static_cast<int>(KF.m_ikp2KF.size());
      const int ip = KF.PushCameraPrior(_iKF1, &m_work);
      if (ip == -1 || static_cast<int>(KF.m_ikp2KF.size()) == Nkp) {
        continue;
      }
      for (int jKF = _iKF2; jKF < nKFs; ++jKF) {
        ++m_iKF2cb[jKF + 1];
      }
    }
  }
}

void GlobalBundleAdjustor::SetCameraPriorMotion(const CameraPriorMotion &Zp) {
  const int nKFs = static_cast<int>(m_KFs.size());
  if (m_ZpLM.Valid()) {
    m_ApLM.MakeMinus();
    m_SAcus[m_ZpLM.m_iKF].Increase3(m_ApLM.m_Arr.m_A, m_ApLM.m_Arr.m_b);
    const int im = m_ZpLM.m_iKF - (nKFs - m_CsLM.Size());
    Camera::Factor::Unitary &SAcm = m_SAcmsLM[im].m_Au;
    SAcm.m_Acm.Increase3(m_ApLM.m_Arm);
    SAcm.m_Amm += m_ApLM.m_Amm;
  }
  m_ZpLM = Zp;
  const int im = m_ZpLM.m_iKF - (nKFs - m_CsLM.Size());
  m_ucmsLM[im] |= GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                  GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                  GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
  m_ApLM.MakeZero();
}

void GlobalBundleAdjustor::PushFeatureMeasurementMatchesFirst(const FRM::Frame &F,
                                                              std::vector<int> &iKF2X,
                                                              std::vector<int> &iX2z) {
  int SNx = 0;
  const int NZ = static_cast<int>(F.m_Zs.size());
  iKF2X.assign(m_KFs.size(), -1);
  for (int iZ = 0; iZ < NZ; ++iZ) {
    const int iKF = F.m_Zs[iZ].m_iKF;
    iKF2X[iKF] = SNx;
    SNx += static_cast<int>(m_KFs[iKF].m_xs.size());
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

void GlobalBundleAdjustor::PushFeatureMeasurementMatchesNext(const FRM::Frame &F1,
                                                             const FRM::Frame &F2, const std::vector<int> &iKF2X, const std::vector<int> &iX2z2,
                                                             FRM::MeasurementMatch &Zm) {
  m_izmsTmp.resize(0);
  const int NZ1 = static_cast<int>(F1.m_Zs.size());
  for (int iZ1 = 0; iZ1 < NZ1; ++iZ1) {
    const FRM::Measurement &Z1 = F1.m_Zs[iZ1];
    const int iX = iKF2X[Z1.m_iKF];
    if (iX == -1) {
      continue;
    }
    const int *ix2z2 = iX2z2.data() + iX;
    const int iz11 = Z1.m_iz1, iz12 = Z1.m_iz2;
    for (int iz1 = iz11; iz1 < iz12; ++iz1) {
      const int iz2 = ix2z2[F1.m_zs[iz1].m_ix];
      if (iz2 != -1) {
        m_izmsTmp.push_back(FTR::Measurement::Match(iz1, iz2));
      }
    }
  }
  Zm.PushFeatureMeasurementMatches(m_izmsTmp);
}

int GlobalBundleAdjustor::CountMeasurementsFrame() {
  int SN = 0;
  const int nKFs = int(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    SN += static_cast<int>(m_KFs[iKF].m_Zs.size());
  }
  return SN;
}

int GlobalBundleAdjustor::CountMeasurementsFeature() {
  int SN = 0;
  const int nKFs = int(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    SN += static_cast<int>(m_KFs[iKF].m_zs.size());
  }
  return SN;
}

int GlobalBundleAdjustor::CountMeasurementsPriorCameraPose() {
  int SN = 0;
  return SN;
}

int GlobalBundleAdjustor::CountSchurComplements() {
  return m_SAcus.Size() + CountSchurComplementsOffDiagonal();
}

int GlobalBundleAdjustor::CountSchurComplementsOffDiagonal() {
  return m_iKF2cb.back();
}
