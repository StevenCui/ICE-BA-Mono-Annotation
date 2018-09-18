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
#ifdef CFG_DEBUG_EIGEN
//#define GBA_DEBUG_EIGEN_PCG
#endif
#include "GlobalBundleAdjustor.h"
#include <math.h>  // isfinite

#if defined WIN32 && defined CFG_DEBUG && defined CFG_GROUND_TRUTH
//#define GBA_DEBUG_GROUND_TRUTH_STATE
//#ifdef GBA_DEBUG_GROUND_TRUTH_STATE
//#define GBA_DEBUG_GROUND_TRUTH_STATE_ERROR
//#endif
#endif

#if defined CFG_DEBUG && defined CFG_VERBOSE
//#define GBA_DEBUG_PCG_SAVE_RESIDUAL
//#define GBA_DEBUG_PCG_SAVE_RESULT
//#define GBA_DEBUG_PCG_LOAD_RESULT
#endif

void GlobalBundleAdjustor::UpdateFactors() {
  const float add = UT::Inverse(BA_VARIANCE_REGULARIZATION_DEPTH, BA_WEIGHT_FEATURE);
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) 
  {
    KeyFrame &KF = m_KFs[iKF];
    if (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) 
    {
      KF.m_SAcxzs.MakeZero();
      KF.m_SAps.MakeZero();
      m_SAcus[iKF].MakeZero();
    } 
    else 
    {
      const int NZ = int(KF.m_Zs.size());
      for (int iZ = 0; iZ < NZ; ++iZ) 
      {
        if (m_ucs[KF.m_Zs[iZ].m_iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA)
          KF.m_SAcxzs[iZ].MakeZero();
      }
      const int Np = KF.m_SAps.Size();
      for (int ip = 0; ip < Np; ++ip) 
      {
        if (m_ucs[KF.m_iKFsPrior[ip]] & GBA_FLAG_FRAME_UPDATE_CAMERA) {
          KF.m_SAps[ip].MakeZero();
        }
      }
    }
    if (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_DEPTH) 
    {
      const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
      const int Nx = static_cast<int>(KF.m_xs.size());
      for (int ix = 0; ix < Nx; ++ix) 
      {
        if (uds[ix] & GBA_FLAG_TRACK_UPDATE_DEPTH) 
        {
          KF.m_Axs[ix].MakeZero();
          KF.m_Axs[ix].m_Sadx.m_add.m_a = add;
        }
      }
    }
  }
  const ubyte ucmFlag = GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
  const int Nm = m_CsLM.Size();
  for (int im = 0; im < Nm; ++im) {
    if (m_ucmsLM[im] & ucmFlag) {
      m_SAcmsLM[im].MakeZero();
    }
  }

  UpdateFactorsFeature();

  UpdateFactorsPriorDepth();  
  UpdateFactorsPriorCameraPose();
  UpdateFactorsPriorCameraMotion();

  UpdateFactorsIMU();

  UpdateFactorsFixOrigin();
  UpdateFactorsFixPositionZ();
  UpdateFactorsFixMotion();
}

void GlobalBundleAdjustor::UpdateFactorsFeature() 
{
  //Tr[0] is T_cz_from_cx
  Rigid3D Tr[2];
  //dadx.m_adc is J_Tz^T * J_dj
  FTR::Factor::DDC dadx;

  Camera::Factor::Unitary::CC dAcxx, dAczz;
  Camera::Factor::Binary::CC dAcxz;
  FTR::Factor::Full::U U;
  
  //m_KFs are all the keyframes in GBA
  const int nKFs = int(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) 
  {
    //KF is measurement frame cz
    KeyFrame &KF = m_KFs[iKF];
    //C is T_cz_from_w
    const Rigid3D &C = m_Cs[iKF];
    //if cz is re-linearized last iteration, ucz is false?
    const bool ucz = (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    
    //m_SAcus is triangle matrix in U
    //SAczz contains m_A(J_cz^T * J_cz), m_r(?) and m_b(J_cz^T * e)
    Camera::Factor::Unitary::CC &SAczz = m_SAcus[iKF];
    //m_Zs is the source KeyFrame from which the current measurements are extracted 
    const int NZ = int(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) 
    {
      //Z is the source keyframe measurements
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      //m_iKF is the id of the source camera cx
      const int _iKF = Z.m_iKF;

      //ucx is false if cx was re-linearized last iteration
      //ucr is false if cz or cx was re-linearized last iteration
      const bool ucx = (m_ucs[_iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0, ucr = ucx || ucz;
      if (!ucr && !(m_ucs[_iKF] & GBA_FLAG_FRAME_UPDATE_DEPTH)) {
        continue;
      }
      //m_Cs[_iKF] is T_cx_from_w, Tr is T_cz_from_cx 
      *Tr = C / m_Cs[_iKF];

      //m_iKF2d is the feature id for the source frame cx
      const int id = m_iKF2d[_iKF];
      const Depth::InverseGaussian *_ds = m_ds.data() + id;
      ubyte *_uds = m_uds.data() + id;
      //SAcxx contains m_A(J_cx^T * J_cx), m_r(J_xj^T * e_cx_cz_j) and m_b(J_cx^T * e)
      Camera::Factor::Unitary::CC &SAcxx = m_SAcus[_iKF];
      //SAcxz is J_cx^T * J_cz 
      Camera::Factor::Binary::CC &SAcxz = KF.m_SAcxzs[iZ];

      //_KF is the source frame cx 
      KeyFrame &_KF = m_KFs[_iKF];
      //m_iz1 and m_iz2 is the feature id in measurement frame cz
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) 
      {
        //m_zs is the feature infomation in source frame cx        
        const FTR::Measurement &z = KF.m_zs[iz];
        //ix is the map point id in the source keyframe cx
        const int ix = z.m_ix;
        //ud is false if landmark's depth was re-linearized last iteration
        const bool ud = (_uds[ix] & GBA_FLAG_TRACK_UPDATE_DEPTH) != 0;
        if (!ucr && !ud) {
          continue;
        }
        //m_Azs2 are all Hession blocks that relate to this visual residua, 
        //and include H_d_d, A_cx_cz, A_cx_cx, A_cz_cz, b_cx_e, b_cz_e
        FTR::Factor::Full::A2 &A = KF.m_Azs2[iz];
        if (!ud) {
          dadx = A.m_adx;
        }
        if (!ucx) {
          dAcxx = A.m_Acxx;
        }
        if (!ucr) {
          dAcxz = A.m_Acxz;
        }
        if (!ucz) {
          dAczz = A.m_Aczz;
        }
        //caculate the Jacobian
        FTR::GetFactor<GBA_ME_FUNCTION>(
          BA_WEIGHT_FEATURE, Tr, _KF.m_xs[ix], _ds[ix], 
          C, z, &KF.m_Lzs[iz], &KF.m_Azs1[iz], &A, &U);
        //ud is true that means depth is not changed last iteration
        //m_Axs are Hession blocks that relate to all landmarks in cx
        if (ud) {
          _KF.m_Axs[ix].m_Sadx += A.m_adx;
        } else {
          FTR::Factor::DDC::amb(A.m_adx, dadx, dadx);
          _KF.m_Axs[ix].m_Sadx += dadx;
        }

        if (ucx) {
          SAcxx += A.m_Acxx;
        } else {
          Camera::Factor::Unitary::CC::AmB(A.m_Acxx, dAcxx, dAcxx);
          SAcxx += dAcxx;
        }

        if (ucr) {
          SAcxz += A.m_Acxz;
        } else {
          Camera::Factor::Binary::CC::AmB(A.m_Acxz, dAcxz, dAcxz);
          SAcxz += dAcxz;
        }

        //if ucz is false that means cz has changed last iteration, 
        // let SAczz = 0 and re-caculate all the Jacobian with cz?
        if (ucz) {
          SAczz += A.m_Aczz;
        } else {
          Camera::Factor::Unitary::CC::AmB(A.m_Aczz, dAczz, dAczz);
          SAczz += dAczz;
        }

        m_ucs[_iKF] |= GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
        _uds[ix] |= GBA_FLAG_TRACK_UPDATE_INFORMATION;

      }
    }
  }
}

void GlobalBundleAdjustor::UpdateFactorsPriorCameraPose() {

  Camera::Factor::Unitary::CC Au, dAu;
  Camera::Factor::Binary::CC Ab, dAb;
  CameraPrior::Pose::Factor::Auxiliary U;
  //float dF;
  const int NZ = static_cast<int>(m_Zps.size());

  for (int iZ = 0; iZ < NZ; ++iZ) {
    const CameraPrior::Pose &Z = m_Zps[iZ];
    const bool ucr = (m_ucs[Z.m_iKFr] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    bool uc = ucr;
    const int N = static_cast<int>(Z.m_iKFs.size());
    for (int i = 0; i < N && !uc; ++i) {
      uc = (m_ucs[Z.m_iKFs[i]] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    }
    if (!uc) {
      continue;
    }

    CameraPrior::Pose::Factor &A = m_Aps[iZ];
    A.Swap(m_ApTmp, m_bpTmp);
    m_work.Resize(U.BindSize(N) / sizeof(float));
    U.Bind(m_work.Data(), N);
    if (GBA_PRIOR_CAMERA_POSE_ROBUST) {
      Z.GetError(m_Cs, &A.m_Je.m_e, BA_ANGLE_EPSILON);
      const float F = (Z.GetCost(1.0f, A.m_Je.m_e) - Z.m_xTb) / BA_WEIGHT_FEATURE;
      //const float F = (A.m_F / BA_WEIGHT_PRIOR_CAMERA_POSE - Z.m_xTb) / BA_WEIGHT_FEATURE;
      const float w = BA_WEIGHT_PRIOR_CAMERA_POSE * ME::Weight<GBA_ME_FUNCTION>(F);
      Z.GetFactor(w, m_Cs, &A, &U, BA_ANGLE_EPSILON);
    } else {
      Z.GetFactor(BA_WEIGHT_PRIOR_CAMERA_POSE, m_Cs, &A, &U, BA_ANGLE_EPSILON);
    }

    for (int i = -1, _i = 0; i < N; ++i, ++_i) {
      const int iKF = i == -1 ? Z.m_iKFr : Z.m_iKFs[i];
      Au.Set(A.m_A[_i][_i], A.m_b[_i]);
      uc = (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
      if (uc) {
        m_SAcus[iKF] += Au;
      } else {
        dAu.Set(m_ApTmp[_i][_i], m_bpTmp[_i]);
        Camera::Factor::Unitary::CC::AmB(Au, dAu, dAu);
        m_SAcus[iKF] += dAu;
      }
      if (i == -1) {
        continue;
      }
      A.m_A[0][_i].Get(Ab, Z.m_iKFr > iKF);
      const int iKF1 = std::min(Z.m_iKFr, iKF), iKF2 = std::max(Z.m_iKFr, iKF);
      KeyFrame &KF2 = m_KFs[iKF2];
      const int ip = KF2.SearchPriorKeyFrame(iKF1);
      if (ucr || uc) {
        KF2.m_SAps[ip] += Ab;
      } else {
        m_ApTmp[0][_i].Get(dAb, Z.m_iKFr > iKF);
        Camera::Factor::Binary::CC::AmB(Ab, dAb, dAb);
        KF2.m_SAps[ip] += dAb;
      }
    }
    for (int i1 = 0, _i1 = 1; i1 < N; ++i1, ++_i1) {
      const int iKF1 = Z.m_iKFs[i1];
      const bool uc1 = (m_ucs[iKF1] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
      for (int i2 = i1 + 1, _i2 = i2 + 1; i2 < N; ++i2, ++_i2) {
        const int iKF2 = Z.m_iKFs[i2];

        const bool uc2 = (m_ucs[iKF2] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
        const CameraPrior::Element::CC &_Ab = A.m_A[_i1][_i2];
        KeyFrame &KF2 = m_KFs[iKF2];
        const int ip = KF2.SearchPriorKeyFrame(iKF1);
        if (uc1 || uc2) {
          KF2.m_SAps[ip] += _Ab;
        } else {
          Camera::Factor::Binary::CC::AmB(_Ab, m_ApTmp[_i1][_i2], dAb);
          KF2.m_SAps[ip] += dAb;
        }
      }
    }
  }
}

void GlobalBundleAdjustor::UpdateFactorsPriorCameraMotion() {
  if (m_ZpLM.Invalid()) {
    return;
  }
  const int im = m_ZpLM.m_iKF - m_Cs.Size() + m_CsLM.Size();
  const ubyte ucm = m_ucmsLM[im];
  const bool uc = (ucm & (GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                          GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION)) != 0;
  const bool _ucm = uc || (ucm & (GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                                  GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                                  GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE)) != 0;
  if (!_ucm) {
    return;
  }
  CameraPrior::Motion::Factor::RR dArr;
  CameraPrior::Motion::Factor::RM dArm;
  CameraPrior::Motion::Factor::MM dAmm;
  if (!uc) {
    dArr = m_ApLM.m_Arr;
  }
  if (!ucm) {
    dArm = m_ApLM.m_Arm;
    dAmm = m_ApLM.m_Amm;
  }
  CameraPrior::Motion::Factor::Auxiliary U;
  m_ZpLM.GetFactor(BA_WEIGHT_PRIOR_CAMERA_MOTION, m_CsLM[im], &m_ApLM, &U);
  Camera::Factor::Unitary::CC &SAcc = m_SAcus[m_ZpLM.m_iKF];
  Camera::Factor::Unitary &SAcm = m_SAcmsLM[im].m_Au;
  if (uc) {
    SAcc.Increase3(m_ApLM.m_Arr.m_A, m_ApLM.m_Arr.m_b);
  } else {
    CameraPrior::Motion::Factor::RR::AmB(m_ApLM.m_Arr, dArr, dArr);
    SAcc.Increase3(dArr.m_A, dArr.m_b);
  }
  if (ucm) {
    SAcm.m_Acm.Increase3(m_ApLM.m_Arm);
    SAcm.m_Amm += m_ApLM.m_Amm;
  } else {
    CameraPrior::Motion::Factor::RM::AmB(m_ApLM.m_Arm, dArm, dArm);
    CameraPrior::Motion::Factor::MM::AmB(m_ApLM.m_Amm, dAmm, dAmm);
    SAcm.m_Acm.Increase3(dArm);
    SAcm.m_Amm += dAmm;
  }
}

void GlobalBundleAdjustor::UpdateFactorsPriorDepth() {
  FTR::Factor::DD dadd;

  const int nKFs = int(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_DEPTH)) {
      continue;
    }
    const int id = m_iKF2d[iKF];
    const Depth::InverseGaussian *ds = m_ds.data() + id;
    ubyte *uds = m_uds.data() + id;
    KeyFrame &KF = m_KFs[iKF];
    const Depth::Prior zp(KF.m_d.u(), 1.0f / (BA_VARIANCE_PRIOR_FRAME_DEPTH + KF.m_d.s2()));
    const int Nx = int(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & GBA_FLAG_TRACK_UPDATE_DEPTH)) {
        continue;
      }
      {
        Depth::Prior::Factor &A = KF.m_Apds[ix];
        zp.GetFactor<GBA_ME_FUNCTION>(BA_WEIGHT_PRIOR_DEPTH, ds[ix].u(), A);
        dadd = A;
      }
      KF.m_Axs[ix].m_Sadx.m_add += dadd;
      m_ucs[iKF] |= GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
      uds[ix] |= GBA_FLAG_TRACK_UPDATE_INFORMATION;
    }
  }
}

void GlobalBundleAdjustor::UpdateFactorsIMU() {
  //dAcc1 is m_A(J_c1^T * J_c1), m_r(?) and m_b(J_c1^T * eb12)
  //dAcc2 is m_A(J_c2^T * J_c2), m_r(?) and m_b(J_c2^T * eb12)
  Camera::Factor::Unitary::CC dAcc1, dAcc2;
  //dAcm1 is m_Acm(J_c1^T * J_m1) and m_Amm (J_m1^T * J_m1)
  Camera::Factor::Unitary dAcm1, dAcm2;
  IMU::Delta::Factor::Auxiliary::Global U;

  const ubyte ucFlag = GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                       GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  const ubyte ucmFlag = ucFlag |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;

  const int Nc = int(m_KFs.size());
  //ic1 is the first camera's id, and ic2 is the second camera
  //im1 is the id of the first IMU motion(v,ba,bg), and im2 is the second motion's id
  for (int ic1 = Nc - m_CsLM.Size(), ic2 = ic1 + 1, im1 = 0, im2 = 1; 
       ic2 < Nc;
       ic1 = ic2++, im1 = im2++) 
  {
    //m_us are IMU measurements
    if (m_KFs[ic2].m_us.Empty()) {
      continue;
    }

    //m_ucmsLM is false if the camera pose and motion are re-linearized last iteration?
    const ubyte ucm1 = m_ucmsLM[im1], ucm2 = m_ucmsLM[im2];
    const bool _ucm1 = (ucm1 & ucmFlag) != 0, _ucm2 = (ucm2 & ucmFlag) != 0;
    if (!_ucm1 && !_ucm2) {
      continue;
    }
    //m_AdsLM ?
    IMU::Delta::Factor &A = m_AdsLM[im2];
    //m_SAcmsLM are Hession blocks that relate to all imu motion,
    // and include H_c_m, H_m_m
    Camera::Factor &SAcm1 = m_SAcmsLM[im1], &SAcm2 = m_SAcmsLM[im2];
    const bool uc1 = (ucm1 & ucFlag) != 0, uc2 = (ucm2 & ucFlag) != 0;

    if (!uc1) {
      dAcc1 = A.m_A11.m_Acc;
    }
    if (!_ucm1) {
      dAcm1.m_Acm = A.m_A11.m_Acm;
      dAcm1.m_Amm = A.m_A11.m_Amm;
    }
    if (!uc2) {
      dAcc2 = A.m_A22.m_Acc;
    }
    if (!_ucm2) {
      dAcm2.m_Acm = A.m_A22.m_Acm;
      dAcm2.m_Amm = A.m_A22.m_Amm;
    }
    //caculate the Jacobian of IMU residual
    m_DsLM[im2].GetFactor(BA_WEIGHT_IMU, m_CsLM[im1], m_CsLM[im2], m_K.m_pu,
                          &A, &SAcm2.m_Ab, &U, BA_ANGLE_EPSILON);

    //m_SAcus is triangle matrix in U    
    if (uc1) {
      m_SAcus[ic1] += A.m_A11.m_Acc;
    } else {
      Camera::Factor::Unitary::CC::AmB(A.m_A11.m_Acc, dAcc1, dAcc1);
      m_SAcus[ic1] += dAcc1;
    }

    if (_ucm1) {
      SAcm1.m_Au.m_Acm += A.m_A11.m_Acm;
      SAcm1.m_Au.m_Amm += A.m_A11.m_Amm;
    } else {
      Camera::Factor::Unitary::CM::AmB(A.m_A11.m_Acm, dAcm1.m_Acm, dAcm1.m_Acm);
      Camera::Factor::Unitary::MM::AmB(A.m_A11.m_Amm, dAcm1.m_Amm, dAcm1.m_Amm);
      SAcm1.m_Au += dAcm1;
    }
  
    if (uc2) {
      m_SAcus[ic2] += A.m_A22.m_Acc;
    } else {
      Camera::Factor::Unitary::CC::AmB(A.m_A22.m_Acc, dAcc2, dAcc2);
      m_SAcus[ic2] += dAcc2;
    }

    //dAcm2.m_Acm = J_c2^T * J_m2, size is 6*9
    //dAcm2.m_Amm = J_m2^T * J_m2, size is 9*9
    if (_ucm2) {
      SAcm2.m_Au.m_Acm += A.m_A22.m_Acm;
      SAcm2.m_Au.m_Amm += A.m_A22.m_Amm;
    } else {
      Camera::Factor::Unitary::CM::AmB(A.m_A22.m_Acm, dAcm2.m_Acm, dAcm2.m_Acm);
      Camera::Factor::Unitary::MM::AmB(A.m_A22.m_Amm, dAcm2.m_Amm, dAcm2.m_Amm);
      SAcm2.m_Au += dAcm2;
    }
  }
}

void GlobalBundleAdjustor::UpdateFactorsFixOrigin() {
  const int iKF = 0;
  if (m_KFs[iKF].m_T.m_iFrm != 0 || !(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA)) {
    return;
  }
  //float dF = m_Af->m_F;
  m_Zo.GetFactor(m_Cs[iKF], m_Ao, BA_ANGLE_EPSILON);
  m_SAcus[iKF] += m_Ao.m_A;
}

void GlobalBundleAdjustor::UpdateFactorsFixPositionZ() {
#ifdef CFG_VERBOSE
  int SN = 0;
#endif
  //float dF;
  const Camera::Fix::PositionZ z(BA_WEIGHT_FIX_POSITION_Z, BA_VARIANCE_FIX_POSITION_Z);
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const ubyte uc = m_ucs[iKF];
    if (!(uc & GBA_FLAG_FRAME_UPDATE_CAMERA)) {
      continue;
    }
    Camera::Fix::PositionZ::Factor &A = m_Afps[iKF];
    Camera::Factor::Unitary::CC &SA = m_SAcus[iKF];
    //dF = A.m_F;
    z.GetFactor(m_Cs[iKF].GetPositionZ(), A);
    SA.m_A.m22() += z.m_w;
    SA.m_b.v2() += A.m_b;
    //dF = A.m_F - dF;
    //m_F = dF + m_F;
#ifdef CFG_VERBOSE
    if (m_verbose >= 3)
      ++SN;
#endif
  }
#ifdef CFG_VERBOSE
  if (m_verbose >= 3) {
    UT::Print("  Fix Position Z = %d / %d = %.2f%%\n", SN, nKFs, UT::Percentage(SN, nKFs));
  }
#endif
}

void GlobalBundleAdjustor::UpdateFactorsFixMotion() {
  //float dF;
  const Camera::Fix::Zero zv[3] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_VELOCITY),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_VELOCITY_INITIAL),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_VELOCITY_INVALID)};
  const Camera::Fix::Zero zba[3] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_ACCELERATION),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_ACCELERATION_INITIAL),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_ACCELERATION_INVALID)};
  const Camera::Fix::Zero zbw[3] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_GYROSCOPE),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_GYROSCOPE_INITIAL),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_GYROSCOPE_INVALID)};
  const ubyte ucmFlag = GBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_POSITION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION |
                        GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE;
  const int Nm = m_CsLM.Size();
  for (int im = 0; im < Nm; ++im) {
    const ubyte ucm = m_ucmsLM[im];
    if (!(ucm & ucmFlag)) {
      continue;
    }
    Camera::Fix::Motion::Factor &A = m_AfmsLM[im];
    Camera::Factor::Unitary::MM &SA = m_SAcmsLM[im].m_Au.m_Amm;
    const Camera &C = m_CsLM[im];
    const int i = ucm >> 5;
    if (ucm & GBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY) {
      zv[i].GetFactor(C.m_v, A.m_Av);
    }
    if (ucm & GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION) {
      zba[i].GetFactor(C.m_ba, A.m_Aba);
    }
    if (ucm & GBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE) {
      zbw[i].GetFactor(C.m_bw, A.m_Abw);
    }
    SA.m_A.IncreaseDiagonal012(zv[i].w());   SA.m_b.Increase(0, A.m_Av.m_b);
    SA.m_A.IncreaseDiagonal345(zba[i].w());  SA.m_b.Increase(3, A.m_Aba.m_b);
    SA.m_A.IncreaseDiagonal678(zbw[i].w());  SA.m_b.Increase(6, A.m_Abw.m_b);
  }
}

void GlobalBundleAdjustor::UpdateSchurComplement() {
  //nKFs is the number of all KF in GBA
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) 
  {
    KeyFrame &KF = m_KFs[iKF];
    //m_ucs?
    if (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) 
    {
      //m_SMcus?
      m_SMcus[iKF].MakeZero();
      //m_Zm.m_SMczms?
      KF.m_Zm.m_SMczms.MakeZero();
    } 
    else {
      const int Nk = KF.m_Zm.m_SMczms.Size();
      for (int ik = 0; ik < Nk; ++ik) {
        if (m_ucs[KF.m_iKFsMatch[ik]] & GBA_FLAG_FRAME_UPDATE_CAMERA) {
          KF.m_Zm.m_SMczms[ik].MakeZero();
        }
      }
    }
  }

  int Nd = 0;
  m_idxsTmp1.assign(nKFs, -1);
  int *iKF2X = m_idxsTmp1.data();
  std::vector<int> &iX2d = m_idxsTmp2;
  iX2d.resize(0);
  const float eps = FLT_EPSILON;
  const float epsd = UT::Inverse(BA_VARIANCE_MAX_DEPTH, BA_WEIGHT_FEATURE, eps);

  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION)) {
      continue;
    }
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    const KeyFrame &KF = m_KFs[iKF];
    const int iX = static_cast<int>(iX2d.size()), Nx = static_cast<int>(KF.m_xs.size());
    iKF2X[iKF] = iX;
    iX2d.resize(iX + Nx, -1);
    int *ix2d = iX2d.data() + iX;
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION)) {
        continue;
      } else if (KF.m_Axs[ix].m_Sadx.m_add.m_a > epsd) {
        ix2d[ix] = Nd++;
      } else if (!(uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
        ix2d[ix] = -2;
      }
    }
  }

  const int NmddC = SIMD_FLOAT_CEIL(Nd);
  m_work.Resize(NmddC + Nd * sizeof(xp128f) / sizeof(float));
  float *mdds = m_work.Data();
  xp128f *_mdds = (xp128f *) (mdds + NmddC);

  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }

    const int *ix2d = iX2d.data() + iX;
    const KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      const int id = ix2d[ix];
      if (id >= 0) {
        mdds[id] = KF.m_Axs[ix].m_Sadx.m_add.m_a;
      }
    }
  }

  //SIMD::Add(Nd, UT::Inverse(BA_VARIANCE_REGULARIZATION_DEPTH, BA_WEIGHT_FEATURE), mdds);
  SIMD::Inverse(Nd, mdds);
  for (int id = 0; id < Nd; ++id) {
    _mdds[id].vdup_all_lane(mdds[id]);
  }

  Camera::Factor::Unitary::CC dMcu;
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    const int *ix2d = iX2d.data() + iX;
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    Camera::Factor::Unitary::CC &SMcxx = m_SMcus[iKF];

    const bool uc = (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;

    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      const int id = ix2d[ix];
      if (id == -1) {
        continue;
      }
      FTR::Factor::Full::Source::M2 &M = KF.m_Mxs2[ix];
      if (id >= 0) {
        if (uc || (uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
          FTR::Factor::Full::Source::Marginalize(_mdds[id], KF.m_Axs[ix].m_Sadx, &KF.m_Mxs1[ix], &M);
          SMcxx += M.m_Mcxx;
        } else {
          dMcu = M.m_Mcxx;
          FTR::Factor::Full::Source::Marginalize(_mdds[id], KF.m_Axs[ix].m_Sadx, &KF.m_Mxs1[ix], &M);
          Camera::Factor::Unitary::CC::AmB(M.m_Mcxx, dMcu, dMcu);
          SMcxx += dMcu;
        }

      } else {

        if (!uc) {
          M.m_Mcxx.GetMinus(dMcu);
          SMcxx += dMcu;
        }
      }
    }
  }

  Camera::Factor::Binary::CC dMcb;
  std::vector<int> &iz2adcz = m_idxsTmp3;
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    KeyFrame &KF = m_KFs[iKF];
    Camera::Factor::Unitary::CC &SMczz = m_SMcus[iKF];

    const int Nz = static_cast<int>(KF.m_zs.size());
    m_marksTmp1.assign(Nz, 0);
    iz2adcz.assign(Nz, -1);
    m_adczsTmp.Resize(0);

    const bool ucz = (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    const int NZ = static_cast<int>(KF.m_Zs.size());

    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iKF = Z.m_iKF, _iX = iKF2X[_iKF];
      if (_iX == -1) {
        continue;
      }
      const ubyte *_uds = m_uds.data() + m_iKF2d[_iKF];
      const int *_ix2mdd = iX2d.data() + _iX;
      Camera::Factor::Binary::CC &SMcxz = KF.m_Zm.m_SMczms[Z.m_ik];

      const bool ucx = (m_ucs[_iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0, ucr = ucx || ucz;
      const KeyFrame &_KF = m_KFs[_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = KF.m_zs[iz].m_ix, id = _ix2mdd[ix];
        if (id == -1) {
          continue;
        }
        FTR::Factor::Full::M2 &M = KF.m_Mzs2[iz];
        if (id >= 0) {
          iz2adcz[iz] = m_adczsTmp.Size();
          //adcz = J_d^T * W * J_cz
          LA::ProductVector6f &adcz = m_adczsTmp.Push();
          if (_uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
            FTR::Factor::Full::Marginalize(_mdds[id], _KF.m_Mxs1[ix], KF.m_Azs1[iz], &KF.m_Mzs1[iz], &M, &adcz);
            SMcxz += M.m_Mcxz;
            SMczz += M.m_Mczz;
            m_marksTmp1[iz] = 1;
          } else {
            if (!ucr) {
              dMcb = M.m_Mcxz;
            }
            if (!ucz) {
              dMcu = M.m_Mczz;
            }
            FTR::Factor::Full::Marginalize(_mdds[id], _KF.m_Mxs1[ix], KF.m_Azs1[iz], &KF.m_Mzs1[iz], &M, &adcz);
            if (ucr) {
              SMcxz += M.m_Mcxz;
            } else {
              Camera::Factor::Binary::CC::AmB(M.m_Mcxz, dMcb, dMcb);
              SMcxz += dMcb;
            }
            if (ucz) {
              SMczz += M.m_Mczz;
            } else {
              Camera::Factor::Unitary::CC::AmB(M.m_Mczz, dMcu, dMcu);
              SMczz += dMcu;
            }
          }
        } else {

          if (!ucr) {
            M.m_Mcxz.GetMinus(dMcb);
            SMcxz += dMcb;
          }
          if (!ucz) {
            M.m_Mczz.GetMinus(dMcu);
            SMczz += dMcu;
          }
          iz2adcz[iz] = -2;

        }

      }
    }

    const int Nk = KF.m_Zm.m_SMczms.Size();
    for (int ik = 0; ik < Nk; ++ik) {
      const int _iKF = KF.m_iKFsMatch[ik];

      const KeyFrame &_KF = m_KFs[_iKF];
      Camera::Factor::Binary::CC &SMczm = KF.m_Zm.m_SMczms[ik];

      const bool _ucz = (m_ucs[_iKF] & GBA_FLAG_FRAME_UPDATE_CAMERA) != 0, uczm = _ucz || ucz;
      const int i1 = KF.m_Zm.m_ik2zm[ik], i2 = KF.m_Zm.m_ik2zm[ik + 1];
      for (int i = i1; i < i2; ++i) {
        const FTR::Measurement::Match &izm = KF.m_Zm.m_izms[i];
        const int iadcz = iz2adcz[izm.m_iz2];
        if (iadcz == -1) {
          continue;
        }
        Camera::Factor::Binary::CC &Mczm = KF.m_Zm.m_Mczms[i];
        if (iadcz >= 0) {
          if (uczm || m_marksTmp1[izm.m_iz2]) {
            FTR::Factor::Full::Marginalize(_KF.m_Mzs1[izm.m_iz1], m_adczsTmp[iadcz], Mczm);
            SMczm += Mczm;
          } else {
            dMcb = Mczm;
            FTR::Factor::Full::Marginalize(_KF.m_Mzs1[izm.m_iz1], m_adczsTmp[iadcz], Mczm);
            Camera::Factor::Binary::CC::AmB(Mczm, dMcb, dMcb);
            SMczm += dMcb;
          }

        } else {
          if (!uczm) {
            Mczm.GetMinus(dMcb);
            SMczm += dMcb;
          }

        }

      }
    }

  }
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    const int *ix2d = iX2d.data() + iX;
    const int Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (ix2d[ix] >= 0) {
        uds[ix] &= ~GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      } else if (ix2d[ix] == -2) {
        uds[ix] |= GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      }
    }
  }

}


bool GlobalBundleAdjustor::SolveSchurComplement() {

  const bool scc = SolveSchurComplementPCG();
  if (GBA_EMBEDDED_MOTION_ITERATION) {
    EmbeddedMotionIteration();
  }

  if (!scc) {
    return false;
  }
  return true;
}

bool GlobalBundleAdjustor::SolveSchurComplementPCG() {
  Camera::Factor::Unitary::CC Acc;
  const int pc = 6, pm = 9;
  const int Nc = m_Cs.Size(), Nm = m_CsLM.Size(), Ncp = Nc * pc, Nmp = Nm * pm, N = Ncp + Nmp;
  m_Acus.Resize(Nc);
  m_bs.Resize(N);
  const float ar = UT::Inverse(BA_VARIANCE_REGULARIZATION_ROTATION, BA_WEIGHT_FEATURE);
  const float ap = UT::Inverse(BA_VARIANCE_REGULARIZATION_POSITION, BA_WEIGHT_FEATURE);
  const float av = UT::Inverse(BA_VARIANCE_REGULARIZATION_VELOCITY, BA_WEIGHT_FEATURE);
  const float aba = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_ACCELERATION, BA_WEIGHT_FEATURE);
  const float abw = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE);
  float *b = m_bs.Data();

  for (int ic = 0; ic < Nc; ++ic, b += pc) {
    Camera::Factor::Unitary::CC::AmB(m_SAcus[ic], m_SMcus[ic], Acc);
    Acc.m_A.IncreaseDiagonal(ap, ar);
    Acc.m_A.GetAlignedMatrix6x6f(m_Acus[ic]);
    Acc.m_b.Get(b);
  }
  m_AmusLM.Resize(Nm);

  for (int im = 0; im < Nm; ++im, b += pm) {
    const Camera::Factor::Unitary::MM &Amm = m_SAcmsLM[im].m_Au.m_Amm;
    m_AmusLM[im].Set(Amm.m_A);
    m_AmusLM[im].IncreaseDiagonal(av, aba, abw);
    Amm.m_b.Get(b);
  }
  const int NKp = CountSchurComplementsOffDiagonal();
  m_Acbs.Resize(NKp);
  m_Acbs.MakeZero();

  for (int ic = 0, im = Nm - Nc, iKp = 0; ic < Nc; ++ic, ++im) {
    LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + m_iKF2cb[ic];
    const KeyFrame &KF = m_KFs[ic];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      Acbs[KF.m_Zs[iZ].m_ik] = KF.m_SAcxzs[iZ];
    }
    const int Np = static_cast<int>(KF.m_iKFsPrior.size());
    for (int ip = 0; ip < Np; ++ip) {
      Acbs[KF.m_iAp2kp[ip]] += KF.m_SAps[ip];
    }
    const int Nk = KF.m_Zm.m_SMczms.Size();
    for (int ik = 0; ik < Nk; ++ik) {
      Acbs[ik] -= KF.m_Zm.m_SMczms[ik];
    }
    //if (im >= 1 && !(m_ucmsLM[im] & GBA_FLAG_CAMERA_MOTION_INVALID) &&
    //               !(m_ucmsLM[im - 1] & GBA_FLAG_CAMERA_MOTION_INVALID)) {
    if (im >= 1 && !KF.m_us.Empty()) {
      Acbs[Nk - 1] += m_SAcmsLM[im].m_Ab.m_Acc;
    }
  }
  m_AcbTs.Resize(NKp);
  for (int ikp = 0; ikp < NKp; ++ikp) {
    m_Acbs[ikp].GetTranspose(m_AcbTs[ikp]);
  }
  PrepareConditioner();

  m_xs.Resize(N);
  m_rs.Resize(N);
  m_ps.Resize(N);
  m_zs.Resize(N);
  m_drs.Resize(N);
  m_dxs.Resize(N);

  bool scc = true;
  PCG_TYPE Se2, Se2Pre, e2Max, alpha, beta, F, FMin;
  PCG_TYPE Se2Min, e2MaxMin;
  m_rs = m_bs;

#ifdef CFG_INCREMENTAL_PCG
  m_xsGN.Resize(0);
  m_xsGN.Push((float *) m_xcs.Data(), Ncp);
  m_xsGN.Push((float *) m_xmsLM.Data(), Nmp);
  m_xsGN.MakeMinus();
  m_xs = m_xsGN;
  ApplyA(m_xs, &m_drs);
  m_rs -= m_drs;
#else
  m_xsGN.MakeZero();
  m_xs.MakeZero();
#endif
  ApplyM(m_rs, &m_ps);

  ConvertCameraMotionResiduals(m_rs, m_ps, &Se2, &e2Max);

  ApplyA(m_ps, &m_drs);
  alpha = Se2 / m_ps.Dot(m_drs);

#ifdef _MSC_VER
  if (_finite(alpha)) {
#else
  if (std::isfinite(alpha)) {
#endif  // _MSC_VER
    const PCG_TYPE e2MaxConv[2] = {ME::ChiSquareDistance<3>(BA_PCG_MIN_CONVERGE_PROBABILITY,
                                                            BA_WEIGHT_FEATURE)/* * s*/,
                                   ME::ChiSquareDistance<3>(BA_PCG_MAX_CONVERGE_PROBABILITY,
                                                            BA_WEIGHT_FEATURE)/* * s*/};
    //////////////////////////////////////////////////////////////////////////
    Se2Min = Se2;
    e2MaxMin = e2Max;
    //////////////////////////////////////////////////////////////////////////

    m_drs *= alpha;
    m_rs -= m_drs;
#ifdef CFG_INCREMENTAL_PCG
    m_ps.GetScaled(alpha, m_dxs);
    m_xs += m_dxs;
#else
    m_ps.GetScaled(alpha, m_xs);
#endif
#ifdef CFG_INCREMENTAL_PCG
    Se2Pre = Se2;
    const PCG_TYPE Se2ConvMin = Se2Pre * BA_PCG_MIN_CONVERGE_RESIDUAL_RATIO;
    const PCG_TYPE Se2ConvMax = Se2Pre * BA_PCG_MAX_CONVERGE_RESIDUAL_RATIO;
#else
    const PCG_TYPE Se2ConvMin = Se2 * BA_PCG_MIN_CONVERGE_RESIDUAL_RATIO;
    const PCG_TYPE Se2ConvMax = Se2 * BA_PCG_MAX_CONVERGE_RESIDUAL_RATIO;
#endif
    const int nIters = std::min(N, BA_PCG_MAX_ITERATIONS);
    for (m_iIterPCG = 0; m_iIterPCG < nIters; ++m_iIterPCG) {
      ApplyM(m_rs, &m_zs);
      Se2Pre = Se2;
      ConvertCameraMotionResiduals(m_rs, m_zs, &Se2, &e2Max);

      //////////////////////////////////////////////////////////////////////////
      if (Se2 < Se2Min) {
        Se2Min = Se2;
        e2MaxMin = e2Max;
        m_xsGN = m_xs;
      }
      //////////////////////////////////////////////////////////////////////////

      beta = Se2 / Se2Pre;


      const int i = (Se2Min <= Se2ConvMin && m_iIterPCG >= BA_PCG_MIN_ITERATIONS) ? 0 : 1;
      if (Se2 == 0.0f || e2MaxMin < e2MaxConv[i]) {
        scc = true;
        break;
      } else if (Se2 > Se2ConvMax) {
        scc = false;
        break;
      }

      m_ps *= beta;
      m_ps += m_zs;
      ApplyA(m_ps, &m_drs);
      alpha = Se2 / m_ps.Dot(m_drs);
#ifdef _MSC_VER
      if (!_finite(alpha)) {
#else
      if (!std::isfinite(alpha)) {
#endif  // _MSC_VER
        scc = false;
        break;
      }
      m_drs *= alpha;
      m_rs -= m_drs;
      m_ps.GetScaled(alpha, m_dxs);
      m_xs += m_dxs;
    }
  } else {
    m_iIterPCG = 0;
  }

  m_xsGN.MakeMinus();
  ConvertCameraUpdates(m_xsGN.Data(), &m_xp2s, &m_xr2s);
  ConvertMotionUpdates(m_xsGN.Data() + Ncp, &m_xv2s, &m_xba2s, &m_xbw2s);

  return scc;
}



void GlobalBundleAdjustor::PrepareConditioner() {
  const int pc = 6, pm = 9;
  //const float eps = 0.0f;
  const float eps = FLT_EPSILON;
  const float epsr = UT::Inverse(BA_VARIANCE_MAX_ROTATION, BA_WEIGHT_FEATURE, eps);
  const float epsp = UT::Inverse(BA_VARIANCE_MAX_POSITION, BA_WEIGHT_FEATURE, eps);
  const float epsv = UT::Inverse(BA_VARIANCE_MAX_VELOCITY, BA_WEIGHT_FEATURE, eps);
  const float epsba = UT::Inverse(BA_VARIANCE_MAX_BIAS_ACCELERATION, BA_WEIGHT_FEATURE, eps);
  const float epsbw = UT::Inverse(BA_VARIANCE_MAX_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE, eps);
  const float epsc[pc] = {epsp, epsp, epsp, epsr, epsr, epsr};
  const float epsm[pm] = {epsv, epsv, epsv, epsba, epsba, epsba, epsbw, epsbw, epsbw};
  const int Nb = GBA_PCG_CONDITIONER_BAND, Nc = m_Cs.Size(), Nm = m_CsLM.Size();
  if (Nb <= 1) {
    m_Mcs.Resize(Nc);
    for (int ic = 0; ic < Nc; ++ic) {

      m_Mcs[ic].Set(m_SAcus[ic], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsc);
      //m_Mcs[ic].Set(m_Acus[ic], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsc);
    }
    m_MmsLM.Resize(Nm);
    for (int im = 0; im < Nm; ++im) {
      m_MmsLM[im].Set(m_AmusLM[im], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsm);
    }
    return;
  }

  m_Mcc.Resize(Nc, Nb);       m_MccT.Resize(Nc, Nb);
  m_McmLM.Resize(Nm, 2);      m_McmTLM.Resize(Nm, 2);
  m_MmcLM.Resize(Nm, Nb - 1); m_MmcTLM.Resize(Nm, Nb - 1);
  m_MmmLM.Resize(Nm, 2);      m_MmmTLM.Resize(Nm, 2);
  m_Mcc.MakeZero();
  for (int ic = 0; ic < Nc; ++ic) {
    m_Mcc[ic][0] = m_Acus[ic];
    //sc.Set(scs[ic]);
    //m_Acus[ic].GetScaledColumn(sc, m_Mcc[ic][0]);
    //m_Mcc[ic][0].ScaleRow(scs[ic]);
    const LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + m_iKF2cb[ic];
    const KeyFrame &KF = m_KFs[ic];
    const int Nkp = static_cast<int>(KF.m_ikp2KF.size());
    for (int ikp = 0; ikp < Nkp; ++ikp) {
      const int _ic = KF.m_ikp2KF[ikp], ib = ic - _ic;
      if (ib >= Nb) {
        continue;
      }
      m_Mcc[_ic][ib] = Acbs[ikp];
      //Acbs[ikp].GetScaledColumn(sc, m_Mcc[_ic][ib]);
      //m_Mcc[_ic][ib].ScaleRow(scs[_ic]);
    }
  }
  for (int ic = Nc - Nm, im = 0; ic < Nc; ++ic, ++im) {
    m_McmLM[im][0] = m_SAcmsLM[im].m_Au.m_Acm;
    m_MmmLM[im][0] = m_AmusLM[im];
    //if (im == 0) {
    //  sm.Set(sms[im]);
    //}
    //m_SAcmsLM[im].m_Au.m_Acm.GetScaledColumn(sm, m_McmLM[im][0]);
    //m_McmLM[im][0].ScaleRow(scs[ic]);
    //m_AmusLM[im].GetScaledColumn(sm, m_MmmLM[im][0]);
    //m_MmmLM[im][0].ScaleRow(sms[im]);
    const int _im = im + 1;
    if (_im == Nm) {
      continue;
    }
    const Camera::Factor::Binary &Ab = m_SAcmsLM[_im].m_Ab;
    m_McmLM[im][1] = Ab.m_Acm;
    m_MmmLM[im][1] = Ab.m_Amm;
    //sm.Set(sms[_im]);
    //Ab.m_Acm.GetScaledColumn(sm, m_McmLM[im][1]);
    //m_McmLM[im][1].ScaleRow(scs[ic]);
    //Ab.m_Amm.GetScaledColumn(sm, m_MmmLM[im][1]);
    //m_MmmLM[im][1].ScaleRow(sms[im]);
    
    LA::AlignedMatrix9x6f *Amcs = m_MmcLM[im];
    Amcs[0] = Ab.m_Amc;
    //Ab.m_Amc.GetScaledColumn(scs[ic + 1], Amcs[0]);
    //Amcs[0].ScaleRow(sms[im]);
    const int Nbm = (im + Nb > Nm ? Nm - im : Nb) - 1;
    for (int ib = 1; ib < Nbm; ++ib) {
      Amcs[ib].MakeZero();
    }
  }

  AlignedVector<LA::AlignedMatrix6x6f> AccsT;
  AlignedVector<LA::AlignedMatrix9x6f> AcmsTLM;
  AlignedVector<LA::AlignedMatrix6x9f> AmcsTLM;
  AlignedVector<LA::AlignedMatrix9x9f> AmmsTLM;
  m_work.Resize((AccsT.BindSize(Nb) + AcmsTLM.BindSize(2) +
                 AmcsTLM.BindSize(Nb) + AmmsTLM.BindSize(2)) / sizeof(float));
  AccsT.Bind(m_work.Data(), Nb);
  AcmsTLM.Bind(AccsT.BindNext(), 2);
  AmcsTLM.Bind(AcmsTLM.BindNext(), Nb);
  AmmsTLM.Bind(AmcsTLM.BindNext(), Nb);
  for (int ic = 0, im = Nm - Nc; ic < Nc; ++ic, ++im) {
    LA::AlignedMatrix6x6f *Mccs = m_Mcc[ic], *MccsT = m_MccT[ic];
    LA::AlignedMatrix6x9f *McmsLM = im >= 0 ? m_McmLM[im] : NULL;
    LA::AlignedMatrix9x6f *McmsTLM = im >= 0 ? m_McmTLM[im] : NULL;
    LA::AlignedMatrix9x6f *MmcsLM = im >= 0 ? m_MmcLM[im] : NULL;
    LA::AlignedMatrix6x9f *MmcsTLM = im >= 0 ? m_MmcTLM[im] : NULL;
    LA::AlignedMatrix9x9f *MmmsLM = im >= 0 ? m_MmmLM[im] : NULL;
    LA::AlignedMatrix9x9f *MmmsTLM = im >= 0 ? m_MmmTLM[im] : NULL;
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = im >= 0 ? (ic + 1 == Nc ? 1 : 2) : 0;
    const int Nbmc = im >= 0 ? Nbcc - 1 : 0;
    const int Nbmm = Nbcm;

    LA::AlignedMatrix6x6f &Mcc = Mccs[0];
    if (Mcc.InverseLDL(epsc)) {
      MccsT[0] = Mcc;
      Mcc.MakeMinus();
      for (int ib = 1; ib < Nbcc; ++ib) {
        Mccs[ib].GetTranspose(AccsT[ib]);
        LA::AlignedMatrix6x6f::ABT(Mcc, AccsT[ib], Mccs[ib]);
        Mccs[ib].GetTranspose(MccsT[ib]);
      }
      for (int ib = 0; ib < Nbcm; ++ib) {
        McmsLM[ib].GetTranspose(AcmsTLM[ib]);
        LA::AlignedMatrix9x6f::ABT(Mcc, AcmsTLM[ib], McmsLM[ib]);
        McmsLM[ib].GetTranspose(McmsTLM[ib]);
      }
      for (int ib = 1; ib < Nbcc; ++ib) {
        const LA::AlignedMatrix6x6f &MccT = MccsT[ib];
        const int _ic = ic + ib;
        LA::AlignedMatrix6x6f *_Mccs = m_Mcc[_ic] - ib;
        LA::AlignedMatrix6x6f::AddABTToUpper(MccT, AccsT[ib], _Mccs[ib]);

        for (int jb = ib + 1; jb < Nbcc; ++jb) {
          LA::AlignedMatrix6x6f::AddABTTo(MccT, AccsT[jb], _Mccs[jb]);
        }
        if (ib == 1 && im >= 0) {
          LA::AlignedMatrix9x6f::AddABTTo(MccT, AcmsTLM[ib], m_McmLM[im + ib][0]);
        }
      }
      for (int ib = 0; ib < Nbcm; ++ib) {
        const LA::AlignedMatrix9x6f &McmT = McmsTLM[ib];
        const int _im = im + ib;
        const LA::AlignedMatrix6x6f *_AccsT = AccsT.Data() + 1;
        LA::AlignedMatrix9x6f *_MmcsLM = m_MmcLM[_im] - ib;
        for (int jb = ib; jb < Nbmc; ++jb) {
          LA::AlignedMatrix9x6f::AddABTTo(McmT, _AccsT[jb], _MmcsLM[jb]);
        }
        LA::AlignedMatrix9x9f *_MmmsLM = m_MmmLM[_im];

        LA::AlignedMatrix9x9f::AddABTToUpper(McmT, AcmsTLM[ib], _MmmsLM[0]);
        if (ib == 0 && Nbcm == 2) {
          LA::AlignedMatrix9x9f::AddABTTo(McmT, AcmsTLM[1], _MmmsLM[1]);
        }
      }
    } else {
      for (int ib = 0; ib < Nbcc; ++ib) {
        Mccs[ib].MakeZero();
        MccsT[ib].MakeZero();
      }
      for (int ib = 0; ib < Nbcm; ++ib) {
        McmsLM[ib].MakeZero();
        McmsTLM[ib].MakeZero();
      }
    }
    if (im < 0) {
      continue;
    }
    LA::AlignedMatrix9x9f &Mmm = MmmsLM[0];

    if (Mmm.InverseLDL(epsm)) {
      MmmsTLM[0] = Mmm;
      Mmm.MakeMinus();

      for (int ib = 0; ib < Nbmc; ++ib) {
        MmcsLM[ib].GetTranspose(AmcsTLM[ib]);
        LA::AlignedMatrix9x9f::ABT(Mmm, AmcsTLM[ib], MmcsLM[ib]);
        MmcsLM[ib].GetTranspose(MmcsTLM[ib]);
      }
      if (Nbmm == 2) {
        MmmsLM[1].GetTranspose(AmmsTLM[1]);
        LA::AlignedMatrix9x9f::ABT(Mmm, AmmsTLM[1], MmmsLM[1]);
        MmmsLM[1].GetTranspose(MmmsTLM[1]);
      }
      for (int ib = 0; ib < Nbmc; ++ib) {
        const LA::AlignedMatrix6x9f &MmcT = MmcsTLM[ib];
        const int _ic = ic + ib + 1;
        LA::AlignedMatrix6x6f *_Mccs = m_Mcc[_ic] - ib;
        LA::AlignedMatrix6x9f::AddABTToUpper(MmcT, AmcsTLM[ib], _Mccs[ib]);

        for (int jb = ib + 1; jb < Nbmc; ++jb) {
          LA::AlignedMatrix6x9f::AddABTTo(MmcT, AmcsTLM[jb], _Mccs[jb]);
        }
        if (ib == 0 && Nbmm == 2) {
          LA::AlignedMatrix9x9f::AddABTTo(MmcT, AmmsTLM[1], m_McmLM[im + ib + 1][0]);
        }
      }
      if (Nbmm == 2) {
        const LA::AlignedMatrix9x9f &MmmT = MmmsTLM[1];
        const int _im = im + 1;
        LA::AlignedMatrix9x6f *_Mmcs = m_MmcLM[_im] - 1;
        for (int jb = 1; jb < Nbmc; ++jb) {
          LA::AlignedMatrix9x9f::AddABTTo(MmmT, AmcsTLM[jb], _Mmcs[jb]);
        }
        LA::AlignedMatrix9x9f::AddABTToUpper(MmmT, AmmsTLM[1], m_MmmLM[_im][0]);
      }
    } else {
      for (int ib = 0; ib < Nbmc; ++ib) {
        MmcsLM[ib].MakeZero();
        MmcsTLM[ib].MakeZero();
      }
      for (int ib = 0; ib < Nbmm; ++ib) {
        MmmsLM[ib].MakeZero();
        MmmsTLM[ib].MakeZero();
      }
    }
  }
}

void GlobalBundleAdjustor::ApplyM(const LA::AlignedVectorX<PCG_TYPE> &xs,
                                  LA::AlignedVectorX<PCG_TYPE> *Mxs) {
  const int Nb = GBA_PCG_CONDITIONER_BAND, Nc = m_Cs.Size(), Nm = m_CsLM.Size();
  if (Nb <= 1) {
#ifdef CFG_PCG_FULL_BLOCK
    LA::ProductVector6f xc;
    LA::AlignedVector9f xm;
#else
    LA::AlignedVector3f xp, xr, xv, xba, xbw;
#endif
    Mxs->Resize(xs.Size());
    const LA::Vector6<PCG_TYPE> *xcs = (LA::Vector6<PCG_TYPE> *) xs.Data();
    LA::Vector6<PCG_TYPE> *Mxcs = (LA::Vector6<PCG_TYPE> *) Mxs->Data();
    for (int ic = 0; ic < Nc; ++ic) {
#ifdef CFG_PCG_FULL_BLOCK
      xc.Set(xcs[ic]);
      m_Mcs[ic].Apply(xc, (PCG_TYPE *) &Mxcs[ic]);
#else
      xcs[ic].Get(xp, xr);
      m_Mcs[ic].Apply(xp, xr, (PCG_TYPE *) &Mxcs[ic]);
#endif
    }
    const LA::Vector9<PCG_TYPE> *xms = (LA::Vector9<PCG_TYPE> *) (xcs + Nc);
    LA::Vector9<PCG_TYPE> *Mxms = (LA::Vector9<PCG_TYPE> *) (Mxcs + Nc);
    const int Nm = m_CsLM.Size();
    for (int im = 0; im < Nm; ++im) {
#ifdef CFG_PCG_FULL_BLOCK
      xm.Set(xms[im]);
      m_MmsLM[im].Apply(xm, (PCG_TYPE *) &Mxms[im]);
#else
      xms[im].Get(xv, xba, xbw);
      m_MmsLM[im].Apply(xv, xba, xbw, (PCG_TYPE *) &Mxms[im]);
#endif
    }
    return;
  }
  
  LA::ProductVector6f bc;
  LA::AlignedVector9f bm;
  Mxs->Set(xs);
  //xs.GetScaled(m_ss, *Mxs);
  LA::Vector6f *bcs = (LA::Vector6f *) Mxs->Data();
  LA::Vector9f *bmsLM = (LA::Vector9f *) (bcs + Nc);

  for (int ic = 0, im = Nm - Nc; ic < Nc; ++ic, ++im) {
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = im >= 0 ? (ic + 1 == Nc ? 1 : 2) : 0;
    const int Nbmc = im >= 0 ? Nbcc - 1 : 0;
    const int Nbmm = Nbcm;
    
    bc.Set(bcs[ic]);
    const LA::AlignedMatrix6x6f *MccsT = m_MccT[ic];
    for (int ib = 1; ib < Nbcc; ++ib) {
      LA::AlignedMatrix6x6f::AddAbTo<float>(MccsT[ib], bc, bcs[ic + ib]);
    }
    if (im >= 0) {
      const LA::AlignedMatrix9x6f *McmsTLM = m_McmTLM[im];
      for (int ib = 0; ib < Nbcm; ++ib) {
        LA::AlignedMatrix9x6f::AddAbTo<float>(McmsTLM[ib], bc, bmsLM[im + ib]);
      }
    }
    LA::AlignedMatrix6x6f::Ab<float>(MccsT[0], bc, bcs[ic]);

    if (im >= 0) {
      bm.Set(bmsLM[im]);
      const LA::AlignedMatrix6x9f *MmcsTLM = m_MmcTLM[im];
      for (int ib = 0; ib < Nbmc; ++ib) {
        LA::AlignedMatrix6x9f::AddAbTo<float>(MmcsTLM[ib], bm, bcs[im + ib + 1]);
      }
      const LA::AlignedMatrix9x9f *MmmsTLM = m_MmmTLM[im];
      if (Nbmm == 2) {
        LA::AlignedMatrix9x9f::AddAbTo<float>(MmmsTLM[1], bm, bmsLM[im + 1]);
      }
      LA::AlignedMatrix9x9f::Ab<float>(MmmsTLM[0], bm, bmsLM[im]);
    }
  }
  m_bcs.Resize(Nc);
  m_bmsLM.Resize(Nm);
  for (int ic = Nc - 1, im = Nm - 1; ic >= 0; --ic, --im) {
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = im >= 0 ? (ic + 1 == Nc ? 1 : 2) : 0;
    const int Nbmc = im >= 0 ? Nbcc - 1 : 0;
    const int Nbmm = Nbcm;

    if (im >= 0) {
      float *_bm = bmsLM[im];
      const LA::AlignedMatrix9x6f *MmcsLM = m_MmcLM[im];
      for (int ib = 0; ib < Nbmc; ++ib) {
        LA::AlignedMatrix9x6f::AddAbTo(MmcsLM[ib], m_bcs[ic + ib + 1], _bm);
      }
      if (Nbmm == 2) {
        LA::AlignedMatrix9x9f::AddAbTo(m_MmmLM[im][1], m_bmsLM[im + 1], _bm);
      }
      m_bmsLM[im].Set(_bm);
    }

    float *_bc = bcs[ic];
    const LA::AlignedMatrix6x6f *Mccs = m_Mcc[ic];
    for (int ib = 1; ib < Nbcc; ++ib) {
      LA::AlignedMatrix6x6f::AddAbTo(Mccs[ib], m_bcs[ic + ib], _bc);
    }
    const LA::AlignedMatrix6x9f *McmsLM = m_McmLM[im];
    for (int ib = 0; ib < Nbcm; ++ib) {
      LA::AlignedMatrix6x9f::AddAbTo(McmsLM[ib], m_bmsLM[im + ib], _bc);
    }
    m_bcs[ic].Set(_bc);
   }
  //Mxs->Scale(m_ss);
}

void GlobalBundleAdjustor::ApplyA(const LA::AlignedVectorX<PCG_TYPE> &xs,
                                  LA::AlignedVectorX<PCG_TYPE> *Axs) {
  const LA::Vector6<PCG_TYPE> *xcs = (LA::Vector6<PCG_TYPE> *) xs.Data();
  ConvertCameraUpdates(xcs, &m_xcsP);
  Axs->Resize(xs.Size());
  LA::Vector6<PCG_TYPE> *Axcs = (LA::Vector6<PCG_TYPE> *) Axs->Data();
  const int Nc = int(m_KFs.size());
  for (int ic = 0; ic < Nc; ++ic) {
    LA::AlignedMatrix6x6f::Ab(m_Acus[ic], m_xcsP[ic], (PCG_TYPE *) &Axcs[ic]);
  }
  for (int ic = 0; ic < Nc; ++ic) {
    const int icb = m_iKF2cb[ic];
    const LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + icb;
    const LA::AlignedMatrix6x6f *AcbTs = m_AcbTs.Data() + icb;
    const LA::ProductVector6f &xc = m_xcsP[ic];
    PCG_TYPE *Axc = Axcs[ic];
    const KeyFrame &KF = m_KFs[ic];
    const int Nkp = static_cast<int>(KF.m_ikp2KF.size());
    for (int ikp = 0; ikp < Nkp; ++ikp) {
      const int _ic = KF.m_ikp2KF[ikp];
      LA::AlignedMatrix6x6f::AddAbTo<PCG_TYPE>(Acbs[ikp], xc, (PCG_TYPE *) &Axcs[_ic]);
      LA::AlignedMatrix6x6f::AddAbTo<PCG_TYPE>(AcbTs[ikp], m_xcsP[_ic], Axc);
    }
  }
  //ConvertCameraUpdates(m_Axcs, Axcs);
  const int Nm = m_CsLM.Size();
  const int ic0 = Nc - Nm;
  ApplyAcm(m_xcsP.Data() + ic0, (LA::Vector9<PCG_TYPE> *) (xcs + Nc), Axcs + ic0,
           (LA::Vector9<PCG_TYPE> *) (Axcs + Nc), false,
           m_AmusLM.Size() == Nm ? m_AmusLM.Data() : NULL);
}

void GlobalBundleAdjustor::ApplyAcm(const LA::ProductVector6f *xcs, const LA::Vector9f *xms,
                                    LA::Vector6f *Axcs, LA::Vector9f *Axms, const bool Acc,
                                    const LA::AlignedMatrix9x9f *Amus) {
  LA::AlignedMatrix6x6f A66;
  LA::AlignedMatrix6x9f A69;
  LA::AlignedMatrix9x6f A96;
  LA::AlignedMatrix9x9f A99;
  LA::AlignedVector9f v9[2];
#ifndef CFG_IMU_FULL_COVARIANCE
  LA::AlignedMatrix3x3f A33;
  LA::AlignedMatrix3x9f A39;
  LA::AlignedVector3f v3;
#endif
  const int Nm = m_CsLM.Size();
  for (int im = 0, r = 0; im < Nm; ++im, r = 1 - r) {
    const LA::ProductVector6f &xc = xcs[im];
    LA::AlignedVector9f &xm = v9[r];
    xm.Set(xms[im]);
    float *Axc = Axcs[im], *Axm = Axms[im];
    const Camera::Factor &SAcm = m_SAcmsLM[im];
    if (Amus) {
      LA::AlignedMatrix9x9f::Ab(Amus[im], xm, Axm);
    } else {
      A99.Set(SAcm.m_Au.m_Amm.m_A);
      LA::AlignedMatrix9x9f::Ab(A99, xm, Axm);
    }
    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Au.m_Acm, xm, Axc);
    SAcm.m_Au.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, xc, Axm);
    if (im == 0) {
      continue;
    }
    const int _im = im - 1;
    const LA::ProductVector6f &_xc = xcs[_im];
    const LA::AlignedVector9f &_xm = v9[1 - r];
    float *_Axc = Axcs[_im], *_Axm = Axms[_im];
    if (Acc) {
      LA::AlignedMatrix6x6f::AddAbTo(SAcm.m_Ab.m_Acc, xc, _Axc);
      SAcm.m_Ab.m_Acc.GetTranspose(A66);
      LA::AlignedMatrix6x6f::AddAbTo(A66, _xc, Axc);
    }
#ifdef CFG_IMU_FULL_COVARIANCE
    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Ab.m_Acm, xm, _Axc);
    SAcm.m_Ab.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, _xc, Axm);
    LA::AlignedMatrix9x6f::AddAbTo(SAcm.m_Ab.m_Amc, xc, _Axm);
    SAcm.m_Ab.m_Amc.GetTranspose(A69);
    LA::AlignedMatrix6x9f::AddAbTo(A69, _xm, Axc);
    LA::AlignedMatrix9x9f::AddAbTo(SAcm.m_Ab.m_Amm, xm, _Axm);
    SAcm.m_Ab.m_Amm.GetTranspose(A99);
    LA::AlignedMatrix9x9f::AddAbTo(A99, _xm, Axm);
#else
    xm.Get012(v3);
    LA::AlignedMatrix3x3f::AddAbTo(SAcm.m_Ab.m_Acm.m_Arv, v3, &_Axc.v3());
    LA::AlignedMatrix9x3f::AddAbTo(SAcm.m_Ab.m_Amm.m_Amv, v3, _Axm);
    SAcm.m_Ab.m_Acm.m_Arv.GetTranspose(A33);
    _xc.Get345(v3);
    LA::AlignedMatrix3x3f::AddAbTo(A33, v3, &Axm.v0());
    LA::AlignedMatrix9x6f::AddAbTo(SAcm.m_Ab.m_Amc, xc, _Axm);
    SAcm.m_Ab.m_Amc.GetTranspose(A69);
    LA::AlignedMatrix6x9f::AddAbTo(A69, _xm, Axc);
    SAcm.m_Ab.m_Amm.m_Amv.GetTranspose(A39);
    LA::AlignedMatrix3x9f::AddAbTo(A39, _xm, &Axm.v0());
    for (int i = 3; i < 9; ++i) {
      const float a = i < 6 ? SAcm.m_Ab.m_Amm.m_Ababa : SAcm.m_Ab.m_Amm.m_Abwbw;
      _Axm[i] = a * xm[i] + _Axm[i];
      Axm[i] = a * _xm[i] + Axm[i];
    }
#endif
  }
}

void GlobalBundleAdjustor::ApplyAcm(const LA::ProductVector6f *xcs, const LA::Vector9d *xms,
                                    LA::Vector6d *Axcs, LA::Vector9d *Axms, const bool Acc,
                                    const LA::AlignedMatrix9x9f *Amus) {
  LA::AlignedMatrix6x6f A66;
  LA::AlignedMatrix6x9f A69;
  LA::AlignedMatrix9x6f A96;
  LA::AlignedMatrix9x9f A99;
  LA::AlignedVector9f v9[2];
#ifndef CFG_IMU_FULL_COVARIANCE
  LA::AlignedMatrix3x3f A33;
  LA::AlignedMatrix3x9f A39;
  LA::AlignedVector3f v3;
#endif
  const int Nm = m_CsLM.Size();
  for (int im = 0, r = 0; im < Nm; ++im, r = 1 - r) {
    const LA::ProductVector6f &xc = xcs[im];
    LA::AlignedVector9f &xm = v9[r];
    xm.Set(xms[im]);
    double *Axc = Axcs[im], *Axm = Axms[im];
    const Camera::Factor &SAcm = m_SAcmsLM[im];
    if (Amus) {
      LA::AlignedMatrix9x9f::Ab(Amus[im], xm, Axm);
    } else {
      A99.Set(SAcm.m_Au.m_Amm.m_A);
      LA::AlignedMatrix9x9f::Ab(A99, xm, Axm);
    }
    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Au.m_Acm, xm, Axc);
    SAcm.m_Au.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, xc, Axm);
    if (im == 0) {
      continue;
    }
    const int _im = im - 1;
    const LA::ProductVector6f &_xc = xcs[_im];
    const LA::AlignedVector9f &_xm = v9[1 - r];
    double *_Axc = Axcs[_im], *_Axm = Axms[_im];
    if (Acc) {
      LA::AlignedMatrix6x6f::AddAbTo(SAcm.m_Ab.m_Acc, xc, _Axc);
      SAcm.m_Ab.m_Acc.GetTranspose(A66);
      LA::AlignedMatrix6x6f::AddAbTo(A66, _xc, Axc);
    }
#ifdef CFG_IMU_FULL_COVARIANCE
    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Ab.m_Acm, xm, _Axc);
    SAcm.m_Ab.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, _xc, Axm);
    LA::AlignedMatrix9x6f::AddAbTo(SAcm.m_Ab.m_Amc, xc, _Axm);
    SAcm.m_Ab.m_Amc.GetTranspose(A69);
    LA::AlignedMatrix6x9f::AddAbTo(A69, _xm, Axc);
    LA::AlignedMatrix9x9f::AddAbTo(SAcm.m_Ab.m_Amm, xm, _Axm);
    SAcm.m_Ab.m_Amm.GetTranspose(A99);
    LA::AlignedMatrix9x9f::AddAbTo(A99, _xm, Axm);
#else
    xm.Get012(v3);
    LA::AlignedMatrix3x3f::AddAbTo(SAcm.m_Ab.m_Acm.m_Arv, v3, &_Axc.v3());
    LA::AlignedMatrix9x3f::AddAbTo(SAcm.m_Ab.m_Amm.m_Amv, v3, _Axm);
    SAcm.m_Ab.m_Acm.m_Arv.GetTranspose(A33);
    _xc.Get345(v3);
    LA::AlignedMatrix3x3f::AddAbTo(A33, v3, &Axm.v0());
    LA::AlignedMatrix9x6f::AddAbTo(SAcm.m_Ab.m_Amc, xc, _Axm);
    SAcm.m_Ab.m_Amc.GetTranspose(A69);
    LA::AlignedMatrix6x9f::AddAbTo(A69, _xm, Axc);
    SAcm.m_Ab.m_Amm.m_Amv.GetTranspose(A39);
    LA::AlignedMatrix3x9f::AddAbTo(A39, _xm, &Axm.v0());
    for (int i = 3; i < 9; ++i) {
      const float a = i < 6 ? SAcm.m_Ab.m_Amm.m_Ababa : SAcm.m_Ab.m_Amm.m_Abwbw;
      _Axm[i] = a * xm[i] + _Axm[i];
      Axm[i] = a * _xm[i] + Axm[i];
    }
#endif
  }
}

GlobalBundleAdjustor::Residual GlobalBundleAdjustor::ComputeResidual(
  const LA::AlignedVectorX<PCG_TYPE> &xs, const bool minus) {
  Residual R;
  const int N = xs.Size();
  m_work.Resize(N);
  LA::AlignedVectorXf rs(m_work.Data(), N, false);
  ApplyA(xs, &rs);
  if (minus) {
    R.m_F = xs.Dot(rs) / 2 - xs.Dot(m_bs);
    rs -= m_bs;
  } else {
    R.m_F = xs.Dot(rs) / 2 + xs.Dot(m_bs);
    rs += m_bs;
  }
  R.m_r2 = rs.SquaredLength();
  return R;
}

void GlobalBundleAdjustor::SolveBackSubstitution() {
  const ubyte ucFlag = GBA_FLAG_FRAME_UPDATE_DELTA | GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const KeyFrame &KF = m_KFs[iKF];
    //const int Nx = static_cast<int>(KF.m_xs.size());
    //if (Nx == 0) {
    //  continue;
    //}
    const bool dc = m_xr2s[iKF] > BA_BACK_SUBSTITUTE_ROTATION ||
                    m_xp2s[iKF] > BA_BACK_SUBSTITUTE_POSITION;
    if (dc) {
      m_ucs[iKF] |= GBA_FLAG_FRAME_UPDATE_DELTA;
    } else {
      m_ucs[iKF] &= ~GBA_FLAG_FRAME_UPDATE_DELTA;
    }
    if (!(m_ucs[iKF] & ucFlag)) {
      continue;
    }
    const int Nx = static_cast<int>(KF.m_xs.size());
    if (Nx > 0) {
      ubyte *uds = m_uds.data() + m_iKF2d[iKF];
      for (int ix = 0; ix < Nx; ++ix) {
        if (!dc && !(uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION)) {
          continue;
        }
        uds[ix] |= GBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION;
        m_ucs[iKF] |= GBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION;
      }
    }
    if (!dc) {
      continue;
    }
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      if (Z.m_iz1 == Z.m_iz2) {
        continue;
      }
      m_ucs[Z.m_iKF] |= GBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION;
      ubyte *_uds = m_uds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        _uds[KF.m_zs[iz].m_ix] |= GBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION;
      }
    }
  }

  int iX = 0;
  m_iKF2X.assign(nKFs, -1);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION)) {
      continue;
    }
    m_iKF2X[iKF] = iX;
    iX += int(m_KFs[iKF].m_xs.size());
  }
  m_xds.Resize(iX);

  LA::AlignedVector6f xc;
  const LA::Vector6f *xcs = (LA::Vector6f *) m_xsGN.Data();
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    //if (!(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION)) {
    //  continue;
    //}
    const ubyte dx = m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_DELTA;
    if (dx) {
      xc.Set(xcs[iKF]);
    }
    const KeyFrame &KF = m_KFs[iKF];
    if (m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION) {
      const LA::AlignedVector6f *_xc = dx ? &xc : NULL;
      const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
      float *xds = m_xds.Data() + m_iKF2X[iKF];
      const int Nx = static_cast<int>(KF.m_xs.size());
      for (int ix = 0; ix < Nx; ++ix) {
        if (!(uds[ix] & GBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION)) {
          continue;
        } else if (uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
          xds[ix] = 0.0f;
        } else {
          xds[ix] = KF.m_Mxs1[ix].BackSubstitute(_xc);
        }
      }
    }
    if (!dx) {
      continue;
    }
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];

      const ubyte *_uds = m_uds.data() + m_iKF2d[Z.m_iKF];
      float *_xds = m_xds.Data() + m_iKF2X[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = KF.m_zs[iz].m_ix;
        if (_uds[ix] & GBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
          continue;
        }
        _xds[ix] = KF.m_Mzs1[iz].BackSubstitute(xc) + _xds[ix];
      }
    }
  }
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = m_iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    const int id = m_iKF2d[iKF];
    const Depth::InverseGaussian *ds = m_ds.data() + id;
    ubyte *uds = m_uds.data() + id;
    float *xds = m_xds.Data() + iX;
    const int Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & GBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION)) {
        continue;
      }
      xds[ix] = -xds[ix];
      //if (Depth::InverseGaussian::Valid(xds[ix] + ds[ix].u())) {
      //  continue;
      //}
      //xds[ix] = 0.0f;
      //uds[ix] &= ~GBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION;
    }
  }
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucs[iKF] & GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION)) {
      continue;
    }
    m_ucs[iKF] &= ~GBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    const int Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      uds[ix] &= ~GBA_FLAG_TRACK_UPDATE_INFORMATION;
    }
  }

  PushDepthUpdates(m_xds, &m_xsGN);
  m_x2GN = m_xsGN.SquaredLength();

}


bool GlobalBundleAdjustor::EmbeddedMotionIteration() {
  const int pc = 6, pm = 9;
  const int Nc = m_Cs.Size(), Nm = m_CsLM.Size();
  const LA::Vector6f *xcs = ((LA::Vector6f *) m_xsGN.Data()) + Nc - Nm;
  LA::Vector9f *xms = (LA::Vector9f *) (xcs + Nm);
  //const float eps = 0.0f;
  const float eps = FLT_EPSILON;

  const int Nmr = Nc * pm, Nmc = pm + pm;
  LA::AlignedMatrixXf A;
  LA::AlignedVectorXf b, ai;
  AlignedVector<LA::AlignedVector18f> x;
  m_work.Resize(A.BindSize(Nmr, Nmc) + b.BindSize(Nmr) + x.BindSize(Nc) + ai.BindSize(Nmc));
  A.Bind(m_work.Data(), Nmr, Nmc);
  b.Bind(A.BindNext(), Nmr);
  x.Bind(b.BindNext(), Nc);
  ai.Bind(x.BindNext(), Nmc);

  LA::AlignedVector6f xc[2];
  LA::AlignedMatrix9x9f Amm;
  LA::AlignedMatrix9x6f Amc;
  for (int im1 = -1, im2 = 0, imp1 = -pm, imp2 = 0, r = 0; im2 < Nc;
       im1 = im2++, imp1 = imp2, imp2 += pm, r = 1 - r) {
    const Camera::Factor &_A = m_SAcmsLM[im2];
    Amm.Set(_A.m_Au.m_Amm.m_A);
    A.SetBlock(imp2, 0, Amm);
    float *_b = b.Data() + imp2;
    _A.m_Au.m_Amm.m_b.Get(_b);
    xc[r].Set(xcs[im2]);
    _A.m_Au.m_Acm.GetTranspose(Amc);
    LA::AlignedMatrix9x6f::AddAbTo(Amc, xc[r], _b);
    if (im2 == 0) {
      continue;
    }
    A.SetBlock(imp1, pm, _A.m_Ab.m_Amm);
    _A.m_Ab.m_Acm.GetTranspose(Amc);
    LA::AlignedMatrix9x6f::AddAbTo(Amc, xc[1 - r], _b);
    LA::AlignedMatrix9x6f::AddAbTo(_A.m_Ab.m_Amc, xc[r], b.Data() + imp1);
  }

  LA::AlignedVector9f _mi;
  const float epsv = UT::Inverse(BA_VARIANCE_MAX_VELOCITY, BA_WEIGHT_FEATURE, eps);
  const float epsba = UT::Inverse(BA_VARIANCE_MAX_BIAS_ACCELERATION, BA_WEIGHT_FEATURE, eps);
  const float epsbw = UT::Inverse(BA_VARIANCE_MAX_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE, eps);
  const float _eps[] = {epsv, epsv, epsv, epsba, epsba, epsba, epsbw, epsbw, epsbw};
  for (int im = 0, imp = 0; im < Nc; ++im) {
    for (int ip = 0; ip < pm; ++ip, ++imp) {
      float *mi = A[imp];
      float &ni = b[imp];
      const float aii = mi[ip];
      if (aii <= _eps[ip]) {
        memset(mi, 0, sizeof(float) * Nmc);
        ni = 0.0f;
        continue;
      }
      const float mii = 1.0f / aii;
      mi[ip] = mii;
      ai.Set(mi, Nmc);
      ai.MakeMinus(ip + 1);
      SIMD::Multiply(ip + 1, Nmc, mii, mi);
      ni *= mii;

      int jmp = imp + 1;
      for (int jp = ip + 1; jp < pm; ++jp, ++jmp) {
        const float aij = ai[jp];
        SIMD::MultiplyAddTo(jp, Nmc, aij, mi, A[jmp]);
        b[jmp] += aij * ni;
      }
      if (im == Nc - 1) {
        continue;
      }
      const float *_ai = ai.Data() + pm;
      _mi.Set(mi + pm);
      for (int jp = 0; jp < pm; ++jp, ++jmp) {
        const float aij = _ai[jp];
        SIMD::MultiplyAddTo(jp, pm, aij, _mi, A[jmp]);
        b[jmp] += aij * ni;
      }
    }
  }

  for (int im = 0, imp = 0; im < Nc; ++im, imp += pm) {
    memcpy(x[im], b.Data() + imp, 36);
  }
  for (int im = Nc - 1, imp = Nmr - 1, r = im & 1; im >= 0; --im, r = 1 - r) {
    const int _im = im + 1;
    const int _Nmc = _im == Nc ? pm : Nmc;
    float *xi = x[im];
    if (_im < Nc) {
      memcpy(xi + pm, x[_im], 36);
    }
    for (int ip = pm - 1; ip >= 0; --ip, --imp) {
      xi[ip] -= SIMD::Dot(ip + 1, _Nmc, A[imp], xi);
    }
  }
  for (int im = 0; im < Nc; ++im) {
    xms[im].Set(x[im]);
  }

  m_xsGN.MakeMinus(pc * Nc);
  ConvertMotionUpdates((float *) xms, &m_xv2s, &m_xba2s, &m_xbw2s);
  return true;
}

void GlobalBundleAdjustor::EmbeddedPointIteration(const AlignedVector<Rigid3D> &Cs,
                                                  const std::vector<ubyte> &ucs,
                                                  const std::vector<ubyte> &uds,
                                                  std::vector<Depth::InverseGaussian> *ds) {
  std::vector<int> &iKF2X = m_idxsTmp1, &iX2d = m_idxsTmp2;
  const int nKFs = static_cast<int>(m_KFs.size());
  iKF2X.assign(nKFs, -1);
  iX2d.resize(0);

  int Nd = 0;
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(ucs[iKF] & GBA_FLAG_FRAME_UPDATE_DEPTH)) {
      continue;
    }
    const ubyte *_uds = uds.data() + m_iKF2d[iKF];
    const int iX = static_cast<int>(iX2d.size()), Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    iKF2X[iKF] = iX;
    iX2d.resize(iX + Nx, -1);
    int *ix2d = iX2d.data() + iX;
    for (int ix = 0; ix < Nx; ++ix) {
      if (_uds[ix] & GBA_FLAG_TRACK_UPDATE_DEPTH) {
        ix2d[ix] = Nd++;
      }
    }
  }

  int Nt = 0;
  m_idxsTmp3.resize(Nd + Nd + 1);
  int *Nzs = m_idxsTmp3.data(), *id2z = Nzs + Nd;
  memset(Nzs, 0, sizeof(int) * Nd);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const KeyFrame &KF = m_KFs[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iX = iKF2X[Z.m_iKF];
      if (_iX == -1) {
        continue;
      }
      bool t = false;
      const int *_ix2d = iX2d.data() + _iX;
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = KF.m_zs[iz];
        const int id = _ix2d[z.m_ix];
        if (id == -1) {
          continue;
        }
        ++Nzs[id];
        t = true;
      }
      if (t) {
        ++Nt;
      }
    }
  }

  m_t12s.Resize(Nt);

  id2z[0] = 0;
  for (int id = 0, iz = 0; id < Nd; ++id) {
    id2z[id + 1] = id2z[id] + Nzs[id];
  }
  m_zds.resize(id2z[Nd]);

  LA::Vector3f Rx;
  //LA::SymmetricMatrix2x2f W;
  Nt = 0;
  memset(Nzs, 0, sizeof(int) * Nd);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const KeyFrame &KF = m_KFs[iKF];
    const Rigid3D &C = Cs[iKF];
    const int NZ = static_cast<int>(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iX = iKF2X[Z.m_iKF];
      if (_iX == -1) {
        continue;
      }
      const int *_ix2d = iX2d.data() + _iX;
      bool found = false;
      for (int iz = Z.m_iz1; iz < Z.m_iz2 && !found; ++iz) {
        found = _ix2d[KF.m_zs[iz].m_ix] != -1;
      }
      if (!found) {
        continue;
      }
      const Rigid3D T = C / Cs[Z.m_iKF];
      LA::AlignedVector3f *t = m_t12s.Data() + Nt++;
      T.GetTranslation(*t);
      const KeyFrame &_KF = m_KFs[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = KF.m_zs[iz];
        const int id = _ix2d[z.m_ix];
        if (id == -1) {
          continue;
        }
        T.ApplyRotation(_KF.m_xs[z.m_ix].m_x, Rx);
        //z.m_W.GetScaled(KF.m_Lzs[iz].m_wx, W);
        const LA::SymmetricMatrix2x2f &W = z.m_W;
        const int i = ++Nzs[id];
        m_zds[i].Set(*t, Rx, z.m_z, z.m_W);
      }
    }
  }

  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    Depth::InverseGaussian *_ds = ds->data() + m_iKF2d[iKF];
    const int *ix2d = iX2d.data() + iX;
    const int Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      const int id = ix2d[ix];
      if (id == -1) {
        continue;
      }
      Depth::InverseGaussian &d = _ds[ix];
      const Depth::InverseGaussian dBkp = d;
      if (!Depth::Triangulate(BA_WEIGHT_FEATURE, Nzs[id], m_zds.data() + id2z[id], &d, &m_work, true)) {
        d = dBkp;
      }
    }
  }
}
