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
#include "LocalBundleAdjustor.h"
#include "Vector12.h"

void LocalBundleAdjustor::UpdateFactors() {
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    if (m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA) {
      m_SAcusLF[iLF].MakeZero();
    }
    if (m_ucmsLF[iLF]) {
      m_SAcmsLF[iLF].MakeZero();
    }
  }

  const ubyte ucFlag = LBA_FLAG_FRAME_PUSH_TRACK | LBA_FLAG_FRAME_UPDATE_DEPTH;
  const ubyte udFlag = LBA_FLAG_TRACK_PUSH | LBA_FLAG_TRACK_UPDATE_DEPTH;
  const float add = UT::Inverse(BA_VARIANCE_REGULARIZATION_DEPTH, BA_WEIGHT_FEATURE);

  const int nKFs = static_cast<int>(m_KFs.size());  
  for (int iKF = 0; iKF < nKFs; ++iKF) 
  {
    if (!(m_ucsKF[iKF] & ucFlag)) {
      continue;
    }
    const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    KeyFrame &KF = m_KFs[iKF];
    //m_xs?
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) 
    {
      if (!(uds[ix] & udFlag)) {
        continue;
      } else if (uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH) {
        KF.m_Axps[ix].MakeZero();
        KF.m_AxpsST[ix].MakeZero();
        KF.m_Axs[ix].MakeZero();
        const int iST1 = KF.m_ix2ST[ix], iST2 = KF.m_ix2ST[ix + 1];
        KF.m_AxsST.MakeZero(iST1, iST2 - iST1);
        KF.m_Axs[ix].m_Sadd.m_a = add;
        for (int iST = iST1; iST < iST2; ++iST) {
          KF.m_AxsST[iST].m_Sadd.m_a = add;
        }
      } else {
        const int iST1 = KF.m_ix2ST[ix], iST2 = KF.m_ix2ST[ix + 1];
        for (int iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
          KF.m_AxsST[iST].m_Sadd.m_a = add;
        }
      }
    }
  }

  UpdateFactorsFeaturePriorDepth();
  UpdateFactorsPriorCameraMotion();

  UpdateFactorsIMU();

  UpdateFactorsFixPositionZ();
  UpdateFactorsFixMotion();
}

void LocalBundleAdjustor::UpdateFactorsFeaturePriorDepth() {
  UpdateFactorsFeatureLF();
  UpdateFactorsFeatureKF();
  UpdateFactorsPriorDepth();
}

void LocalBundleAdjustor::UpdateFactorsFeatureLF() {
  Rigid3D Tr[2];
  FTR::Factor::DD dadd, daddST;
  Camera::Factor::Unitary::CC dAczz;
  FTR::Factor::FixSource::U U;
  const ubyte ucFlag = LBA_FLAG_FRAME_PUSH_TRACK | LBA_FLAG_FRAME_POP_TRACK |
                       LBA_FLAG_FRAME_UPDATE_DEPTH;
  const ubyte udFlag = LBA_FLAG_TRACK_PUSH | LBA_FLAG_TRACK_POP |
                       LBA_FLAG_TRACK_UPDATE_DEPTH;

  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) 
  {
    LocalFrame &LF = m_LFs[iLF];
    const Rigid3D &C = m_CsLF[iLF].m_T;
    const bool ucz = (m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    //SAczz: 6x6f m_A, float m_r, 6x1f m_b
    Camera::Factor::Unitary::CC &SAczz = m_SAcusLF[iLF];

    // NZ is the number of source frame for current frame
    const int NZ = static_cast<int>(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) 
    {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      //iKF is the source frame id
      const int iKF = Z.m_iKF;
      const bool ucx = (m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0, ucr = ucx || ucz;
      if (!ucr && !(m_ucsKF[iKF] & ucFlag)) {
        continue;
      }
      const bool pushFrm = (m_ucsKF[iKF] & LBA_FLAG_FRAME_PUSH_TRACK) != 0;
      *Tr = C / m_CsKF[iKF];
      //id is the map point id for the source frame in all the points
      const int id = m_iKF2d[iKF];
      ubyte *uds = m_uds.data() + id;
      const Depth::InverseGaussian *ds = m_ds.data() + id;
      //KF is the source frame
      KeyFrame &KF = m_KFs[iKF];
      //iz is the feature id in the measurement frame
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) 
      {
        //z is the measurement in source frame        
        const FTR::Measurement &z = LF.m_zs[iz];
        //ix is the feature id in the source frame
        const int ix = z.m_ix;
        if (!ucr && !(uds[ix] & udFlag)) {
          continue;
        }
        const LocalFrame::SlidingTrack &ST = LF.m_STs[iz];
        const int iST0 = KF.m_ix2ST[ix], iST1 = iST0 + ST.m_ist1, iST2 = iST0 + ST.m_ist2;
        const int Nst = iST2 - iST1;
        FTR::Factor::FixSource::A2 &A = LF.m_Azs2[iz];
        FTR::Factor::FixSource::A3 &AST = LF.m_AzsST[iz];
        const bool ud = (uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH) != 0;
        const bool pushST = pushFrm && (KF.m_usST[iST2 - 1] & LBA_FLAG_TRACK_PUSH);
        if (ucr || ud) {
          if (!ud) {
            dadd = A.m_add;
            daddST = AST.m_add;
          }
          if (!ucz) {
            dAczz = A.m_Aczz;
          }
          //dF = A.m_F;
          FTR::GetFactor<LBA_ME_FUNCTION>(BA_WEIGHT_FEATURE, Tr, KF.m_xs[ix], ds[ix], C, z,
                                          &LF.m_Lzs[iz], &LF.m_Azs1[iz], &A, &U);
          AST.Set(A.m_add, LF.m_Azs1[iz].m_adczA);
          if (Nst >= 1) {
            AST *= 1.0f / Nst;
          }
          LF.m_Nsts[iz] = Nst;
          if (ud) {
            KF.m_Axs[ix].m_Sadd += A.m_add;
            for (int iST = iST1; iST < iST2; ++iST) {
              KF.m_AxsST[iST].m_Sadd += AST.m_add;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          } else {
            FTR::Factor::DD::amb(A.m_add, dadd, dadd);
            KF.m_Axs[ix].m_Sadd += dadd;
            if (pushST) {
              int iST;
              for (iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
                KF.m_AxsST[iST].m_Sadd += AST.m_add;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
              FTR::Factor::DD::amb(AST.m_add, daddST, daddST);
              for (; iST >= iST1; --iST) {
                KF.m_AxsST[iST].m_Sadd += daddST;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
            }
            else {
              FTR::Factor::DD::amb(AST.m_add, daddST, daddST);
              for (int iST = iST1; iST < iST2; ++iST) {
                KF.m_AxsST[iST].m_Sadd += daddST;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
            }
          }
          if (ucz) {
            SAczz += A.m_Aczz;
          } else {
            Camera::Factor::Unitary::CC::AmB(A.m_Aczz, dAczz, dAczz);
            SAczz += dAczz;
          }
          //dF = A.m_F - dF;
          //m_F = dF + m_F;
          m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
          uds[ix] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
        } else if (pushST || LF.m_Nsts[iz] != Nst) {
          if (LF.m_Nsts[iz] == Nst) {
            for (int iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
              KF.m_AxsST[iST].m_Sadd += AST.m_add;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          }
          else {
            daddST = AST.m_add;
            AST.Set(A.m_add, LF.m_Azs1[iz].m_adczA);
            if (Nst >= 1) {
              AST *= 1.0f / Nst;
            }
            LF.m_Nsts[iz] = Nst;
            if (pushST) {
              int iST;
              for (iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
                KF.m_AxsST[iST].m_Sadd += AST.m_add;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
              FTR::Factor::DD::amb(AST.m_add, daddST, daddST);
              for (; iST >= iST1; --iST) {
                KF.m_AxsST[iST].m_Sadd += daddST;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
            }
            else {
              FTR::Factor::DD::amb(AST.m_add, daddST, daddST);
              for (int iST = iST1; iST < iST2; ++iST) {
                KF.m_AxsST[iST].m_Sadd += daddST;
                KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              }
            }
          }
          m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
        }
      }
    }
  }
}

void LocalBundleAdjustor::UpdateFactorsFeatureKF() {
  Rigid3D Tr[2];
  FTR::Factor::DD dadd;
  FTR::Factor::Depth::U U;
  //float dF;
  const ubyte ucFlag = LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION |
                       LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION_KF;
  const ubyte udFlag = LBA_FLAG_TRACK_UPDATE_INFORMATION | LBA_FLAG_TRACK_UPDATE_INFORMATION_KF;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    KeyFrame &KF = m_KFs[iKF];
    const Rigid3D &C = m_CsKF[iKF];
    const bool ucz = (m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    const int NZ = int(KF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = KF.m_Zs[iZ];
      const int _iKF = Z.m_iKF;
      const bool ucx = (m_ucsKF[_iKF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0, ucr = ucx || ucz;
      if (!ucr && !(m_ucsKF[_iKF] & LBA_FLAG_FRAME_UPDATE_DEPTH)) {
        continue;
      }
      *Tr = C / m_CsKF[_iKF];
      const int id = m_iKF2d[_iKF];
      ubyte *_uds = m_uds.data() + id;
      const Depth::InverseGaussian *_ds = m_ds.data() + id;
      KeyFrame &_KF = m_KFs[_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = KF.m_zs[iz];
        const int ix = z.m_ix;
        const bool ud = (_uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH) != 0;
        if (!ucr && !ud) {
          continue;
        }
        FTR::Factor::Depth &A = KF.m_Azs[iz];
        if (!ud) {
          dadd = A.m_add;
        }
        //dF = A.m_F;
        FTR::GetFactor<LBA_ME_FUNCTION>(BA_WEIGHT_FEATURE_KEY_FRAME, Tr, _KF.m_xs[ix], _ds[ix],
                                        C, z, &A, &U);
        if (ud) {
          _KF.m_Axps[ix].m_Sadd += A.m_add;
          _KF.m_Axs[ix].m_Sadd += A.m_add;
        } else {
          FTR::Factor::DD::amb(A.m_add, dadd, dadd);
          _KF.m_Axps[ix].m_Sadd += dadd;
          _KF.m_Axs[ix].m_Sadd += dadd;
        }
        //dF = A.m_F - dF;
        //m_F = dF + m_F;
        m_ucsKF[_iKF] |= ucFlag;
        _uds[ix] |= udFlag;
      }
    }
  }
}

void LocalBundleAdjustor::UpdateFactorsPriorDepth() {
  FTR::Factor::DD dadd, daddST;
  const ubyte ucFlag1 = LBA_FLAG_FRAME_PUSH_TRACK | LBA_FLAG_FRAME_POP_TRACK |
                        LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION_KF;
  const ubyte ucFlag2 = ucFlag1 | LBA_FLAG_FRAME_UPDATE_DEPTH;
  const ubyte udFlag1 = LBA_FLAG_TRACK_PUSH | LBA_FLAG_TRACK_POP |
                        LBA_FLAG_TRACK_UPDATE_INFORMATION_KF;
  const ubyte udFlag2 = udFlag1 | LBA_FLAG_TRACK_UPDATE_DEPTH;
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & ucFlag2)) {
      continue;
    }
    const bool pushFrm = (m_ucsKF[iKF] & LBA_FLAG_FRAME_PUSH_TRACK) != 0;
    m_ucsKF[iKF] &= ~ucFlag1;
    const int id = m_iKF2d[iKF];
    const Depth::InverseGaussian *ds = m_ds.data() + id;
    ubyte *uds = m_uds.data() + id;
    KeyFrame &KF = m_KFs[iKF];
    const Depth::Prior zp(KF.m_d.u(), 1.0f / (BA_VARIANCE_PRIOR_FRAME_DEPTH + KF.m_d.s2()));
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & udFlag2)) {
        continue;
      }
      const ubyte ud = uds[ix];
      uds[ix] &= ~udFlag1;
      const int iST1 = KF.m_ix2ST[ix], iST2 = KF.m_ix2ST[ix + 1], Nst = iST2 - iST1;
      const bool _ud = (ud & LBA_FLAG_TRACK_UPDATE_DEPTH) != 0;
      if (_ud) {
        {
          Depth::Prior::Factor &A = KF.m_Apds[ix];
          //dF = a.m_F;
          zp.GetFactor<LBA_ME_FUNCTION>(BA_WEIGHT_PRIOR_DEPTH, ds[ix].u(), A);
          dadd = A;
        }
        KF.m_Axps[ix].m_Sadd += dadd;
        KF.m_Axs[ix].m_Sadd += dadd;
        //dF = a.m_F - dF;
        //m_F = dF + m_F;
        m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
        uds[ix] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
      }
      FTR::Factor::FixSource::Source::A &AST = KF.m_AxpsST[ix];
      if (Nst == 0) {
        AST = KF.m_Axps[ix];
        KF.m_Nsts[ix] = Nst;
        continue;
      }
      const bool pushST = pushFrm && (KF.m_usST[iST2 - 1] & LBA_FLAG_TRACK_PUSH);
      if (_ud || (ud & LBA_FLAG_TRACK_UPDATE_INFORMATION_KF)) {
        if (!_ud) {
          daddST = AST.m_Sadd;
        }
        AST = KF.m_Axps[ix];
        if (Nst > 1) {
          AST *= 1.0f / Nst;
        }
        KF.m_Nsts[ix] = Nst;
        if (_ud) {
          for (int iST = iST1; iST < iST2; ++iST) {
            KF.m_AxsST[iST].m_Sadd += AST.m_Sadd;
            KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
          }
          if (pushST) {
            for (int iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
              KF.m_usST[iST] &= ~LBA_FLAG_TRACK_PUSH;
            }
          }
        } else {
          if (pushST) {
            int iST;
            for (iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
              KF.m_AxsST[iST].m_Sadd += AST.m_Sadd;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              KF.m_usST[iST] &= ~LBA_FLAG_TRACK_PUSH;
            }
            FTR::Factor::DD::amb(AST.m_Sadd, daddST, daddST);
            for (; iST >= iST1; --iST) {
              KF.m_AxsST[iST].m_Sadd += daddST;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          } else {
            FTR::Factor::DD::amb(AST.m_Sadd, daddST, daddST);
            for (int iST = iST1; iST < iST2; ++iST) {
              KF.m_AxsST[iST].m_Sadd += daddST;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          }
        }
      } else if (pushST || KF.m_Nsts[ix] != Nst) {
        if (KF.m_Nsts[ix] == Nst) {
          for (int iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
            KF.m_AxsST[iST].m_Sadd += AST.m_Sadd;
            KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            KF.m_usST[iST] &= ~LBA_FLAG_TRACK_PUSH;
          }
        } else {
          daddST = AST.m_Sadd;
          AST = KF.m_Axps[ix];
          if (Nst > 1) {
            AST *= 1.0f / Nst;
          }
          KF.m_Nsts[ix] = Nst;
          if (pushST) {
            int iST;
            for (iST = iST2 - 1; iST >= iST1 && (KF.m_usST[iST] & LBA_FLAG_TRACK_PUSH); --iST) {
              KF.m_AxsST[iST].m_Sadd += AST.m_Sadd;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
              KF.m_usST[iST] &= ~LBA_FLAG_TRACK_PUSH;
            }
            FTR::Factor::DD::amb(AST.m_Sadd, daddST, daddST);
            for (; iST >= iST1; --iST) {
              KF.m_AxsST[iST].m_Sadd += daddST;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          } else {
            FTR::Factor::DD::amb(AST.m_Sadd, daddST, daddST);
            for (int iST = iST1; iST < iST2; ++iST) {
              KF.m_AxsST[iST].m_Sadd += daddST;
              KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION;
            }
          }
        }
        m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
      }
    }
  }
}

void LocalBundleAdjustor::UpdateFactorsPriorCameraMotion() {
  const int iLF = m_ic2LF.front();
  const ubyte ucm = m_ucmsLF[iLF];
  if (!ucm) {
    return;
  }
  const bool uc = (ucm & (LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                          LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION)) != 0;
  CameraPrior::Motion::Factor::RR dArr;
  CameraPrior::Motion::Factor::RM dArm;
  CameraPrior::Motion::Factor::MM dAmm;
  if (!uc) {
    dArr = m_ApLF.m_Arr;
  }
  if (!ucm) {
    dArm = m_ApLF.m_Arm;
    dAmm = m_ApLF.m_Amm;
  }
  CameraPrior::Motion::Factor::Auxiliary U;
  m_ZpLF.GetFactor(BA_WEIGHT_PRIOR_CAMERA_MOTION, m_CsLF[iLF], &m_ApLF, &U);
  //m_ApLF.Print();
  Camera::Factor::Unitary::CC &SAcc = m_SAcusLF[iLF];
  Camera::Factor::Unitary &SAcm = m_SAcmsLF[iLF].m_Au;
  if (uc) {
    SAcc.Increase3(m_ApLF.m_Arr.m_A, m_ApLF.m_Arr.m_b);
  } else {
    CameraPrior::Motion::Factor::RR::AmB(m_ApLF.m_Arr, dArr, dArr);
    SAcc.Increase3(dArr.m_A, dArr.m_b);
  }
  if (ucm) {
    SAcm.m_Acm.Increase3(m_ApLF.m_Arm);
    SAcm.m_Amm += m_ApLF.m_Amm;
  } else {
    CameraPrior::Motion::Factor::RM::AmB(m_ApLF.m_Arm, dArm, dArm);
    CameraPrior::Motion::Factor::MM::AmB(m_ApLF.m_Amm, dAmm, dAmm);
    SAcm.m_Acm.Increase3(dArm);
    SAcm.m_Amm += dAmm;
  }
}

void LocalBundleAdjustor::UpdateFactorsIMU() {
  Camera::Factor::Unitary::CC dAcc1, dAcc2;
  Camera::Factor::Unitary dAcm1, dAcm2;
  IMU::Delta::Factor::Auxiliary::Global U;
  //float dF;
  const ubyte ucFlag = LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                       LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int ic1 = 0, ic2 = 1; ic2 < nLFs; ic1 = ic2++) {
    const int iLF1 = m_ic2LF[ic1], iLF2 = m_ic2LF[ic2];
    const ubyte ucm1 = m_ucmsLF[iLF1], ucm2 = m_ucmsLF[iLF2];
    if (!ucm1 && !ucm2) {
      continue;
    }
    IMU::Delta::Factor &A = m_AdsLF[iLF2];
    Camera::Factor &SAcm1 = m_SAcmsLF[iLF1], &SAcm2 = m_SAcmsLF[iLF2];
    const bool uc1 = (ucm1 & ucFlag) != 0, uc2 = (ucm2 & ucFlag) != 0;
    if (!uc1) {
      dAcc1 = A.m_A11.m_Acc;
    }
    if (!ucm1) {
      dAcm1.m_Acm = A.m_A11.m_Acm;
      dAcm1.m_Amm = A.m_A11.m_Amm;
    }
    if (!uc2) {
      dAcc2 = A.m_A22.m_Acc;
    }
    if (!ucm2) {
      dAcm2.m_Acm = A.m_A22.m_Acm;
      dAcm2.m_Amm = A.m_A22.m_Amm;
    }
    //dF = A.m_F;
    m_DsLF[iLF2].GetFactor(BA_WEIGHT_IMU, m_CsLF[iLF1], m_CsLF[iLF2], m_K.m_pu,
                           &A, &SAcm2.m_Ab, &U, BA_ANGLE_EPSILON);
    if (uc1) {
      m_SAcusLF[iLF1] += A.m_A11.m_Acc;
    } else {
      Camera::Factor::Unitary::CC::AmB(A.m_A11.m_Acc, dAcc1, dAcc1);
      m_SAcusLF[iLF1] += dAcc1;
    }
    if (ucm1) {
      SAcm1.m_Au.m_Acm += A.m_A11.m_Acm;
      SAcm1.m_Au.m_Amm += A.m_A11.m_Amm;
    } else {
      Camera::Factor::Unitary::CM::AmB(A.m_A11.m_Acm, dAcm1.m_Acm, dAcm1.m_Acm);
      Camera::Factor::Unitary::MM::AmB(A.m_A11.m_Amm, dAcm1.m_Amm, dAcm1.m_Amm);
      SAcm1.m_Au += dAcm1;
    }
    if (uc2) {
      m_SAcusLF[iLF2] += A.m_A22.m_Acc;
    } else {
      Camera::Factor::Unitary::CC::AmB(A.m_A22.m_Acc, dAcc2, dAcc2);
      m_SAcusLF[iLF2] += dAcc2;
    }
    if (ucm2) {
      SAcm2.m_Au.m_Acm += A.m_A22.m_Acm;
      SAcm2.m_Au.m_Amm += A.m_A22.m_Amm;
    } else {
      Camera::Factor::Unitary::CM::AmB(A.m_A22.m_Acm, dAcm2.m_Acm, dAcm2.m_Acm);
      Camera::Factor::Unitary::MM::AmB(A.m_A22.m_Amm, dAcm2.m_Amm, dAcm2.m_Amm);
      SAcm2.m_Au += dAcm2;
    }
  }
}

void LocalBundleAdjustor::UpdateFactorsFixOrigin() {
  const int iLF = m_ic2LF[0];
  if (m_LFs[iLF].m_T.m_iFrm != 0 || !(m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA)) {
    return;
  }
  m_Zo.GetFactor(m_CsLF[iLF].m_T, m_Ao, BA_ANGLE_EPSILON);
  m_SAcusLF[iLF] += m_Ao.m_A;
}

void LocalBundleAdjustor::UpdateFactorsFixPositionZ() {
  //float dF;
  const Camera::Fix::PositionZ z(BA_WEIGHT_FIX_POSITION_Z, BA_VARIANCE_FIX_POSITION_Z);
  const ubyte ucmFlag = LBA_FLAG_CAMERA_MOTION_UPDATE_ROTATION |
                        LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION;
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    const ubyte ucm = m_ucmsLF[iLF];
    if (!(ucm & ucmFlag)) {
      continue;
    }
    Camera::Fix::PositionZ::Factor &A = m_AfpsLF[iLF];
    Camera::Factor::Unitary::CC &SA = m_SAcusLF[iLF];
    if (ucm & LBA_FLAG_CAMERA_MOTION_UPDATE_POSITION) {
      //dF = A.m_F;
      z.GetFactor(m_CsLF[iLF].m_p.z(), A);
    }
    SA.m_A.m22() += z.m_w;
    SA.m_b.v2() += A.m_b;
  }
}

void LocalBundleAdjustor::UpdateFactorsFixMotion() {
  //float dF;
  const Camera::Fix::Zero zv[2] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_VELOCITY),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_VELOCITY_INITIAL)};
  const Camera::Fix::Zero zba[2] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_ACCELERATION),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_ACCELERATION_INITIAL)};
  const Camera::Fix::Zero zbw[2] = {
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_GYROSCOPE),
        Camera::Fix::Zero(BA_WEIGHT_FIX_MOTION, BA_VARIANCE_FIX_BIAS_GYROSCOPE_INITIAL)};
  const int nLFs = int(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    const ubyte ucm = m_ucmsLF[iLF];
    if (!ucm) {
      continue;
    }
    Camera::Fix::Motion::Factor &A = m_AfmsLF[iLF];
    Camera::Factor::Unitary::MM &SA = m_SAcmsLF[iLF].m_Au.m_Amm;
    const Camera &C = m_CsLF[iLF];
    const int i = m_LFs[iLF].m_T.m_iFrm == 0 ? 1 : 0;
    if (ucm & LBA_FLAG_CAMERA_MOTION_UPDATE_VELOCITY) {
      //dF = A.m_Av.m_F;
      zv[i].GetFactor(C.m_v, A.m_Av);
    }
    if (ucm & LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_ACCELERATION) {
      //dF = A.m_Aba.m_F;
      zba[i].GetFactor(C.m_ba, A.m_Aba);
    }
    if (ucm & LBA_FLAG_CAMERA_MOTION_UPDATE_BIAS_GYROSCOPE) {
      zbw[i].GetFactor(C.m_bw, A.m_Abw);
    }
    SA.m_A.IncreaseDiagonal012(zv[i].w());   SA.m_b.Increase(0, A.m_Av.m_b);
    SA.m_A.IncreaseDiagonal345(zba[i].w());  SA.m_b.Increase(3, A.m_Aba.m_b);
    SA.m_A.IncreaseDiagonal678(zbw[i].w());  SA.m_b.Increase(6, A.m_Abw.m_b);
  }
}

void LocalBundleAdjustor::UpdateSchurComplement() {
  const int nLFs = static_cast<int>(m_LFs.size());
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    LocalFrame &LF = m_LFs[iLF];
    if (m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA) {
      m_SMcusLF[iLF].MakeZero();
      LF.m_Zm.m_SMczms.MakeZero();
      LF.m_SmddsST.MakeZero();
      LF.m_Zm.m_SmddsST.MakeZero();
    } else {
      const int Nk = static_cast<int>(LF.m_iLFsMatch.size());
      for (int ik = 0; ik < Nk; ++ik) {
        if (!(m_ucsLF[LF.m_iLFsMatch[ik]] & LBA_FLAG_FRAME_UPDATE_CAMERA)) {
          continue;
        }
        LF.m_Zm.m_SMczms[ik].MakeZero();
        const int izm = LF.m_Zm.m_ik2zm[ik], Nzm = LF.m_Zm.m_ik2zm[ik + 1] - izm;
        LF.m_Zm.m_SmddsST.MakeZero(izm, Nzm);
      }
      const int NZ = int(LF.m_Zs.size());
      for (int iZ = 0; iZ < NZ; ++iZ) {
        const FRM::Measurement &Z = LF.m_Zs[iZ];
        if (!(m_ucsKF[Z.m_iKF] & LBA_FLAG_FRAME_UPDATE_DEPTH)) {
          continue;
        }
        const ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
        for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
          if (uds[LF.m_zs[iz].m_ix] & LBA_FLAG_FRAME_UPDATE_DEPTH) {
            LF.m_SmddsST[iz].MakeZero();
          }
        }
      }
      const int NI = int(LF.m_Zm.m_Is.size());
      for (int iI = 0; iI < NI; ++iI) {
        const MeasurementMatchLF::Index &I = LF.m_Zm.m_Is[iI];
        if (!(m_ucsKF[I.m_iKF] & LBA_FLAG_FRAME_UPDATE_DEPTH) ||
            (m_ucsLF[LF.m_iLFsMatch[I.m_ik]] & LBA_FLAG_FRAME_UPDATE_CAMERA)) {
          continue;
        }
        const ubyte *uds = m_uds.data() + m_iKF2d[I.m_iKF];
        const int i1 = LF.m_Zm.m_iI2zm[iI], i2 = LF.m_Zm.m_iI2zm[iI + 1];
        for (int i = i1; i < i2; ++i) {
          if (uds[LF.m_zs[LF.m_Zm.m_izms[i].m_iz1].m_ix] & LBA_FLAG_FRAME_UPDATE_DEPTH) {
            LF.m_Zm.m_SmddsST[i] = 0.0f;
          }
        }
      }
    }
  }

  int Nd = 0, NdST = 0;
  const int nKFs = int(m_KFs.size());
  m_idxsTmp1.assign(nKFs + nKFs, -1);
  int *iKF2X = m_idxsTmp1.data(), *iKF2XST = m_idxsTmp1.data() + nKFs;
  std::vector<int> &iX2d = m_idxsTmp2, &iXST2dST = m_idxsTmp3;
  iX2d.resize(0);
  iXST2dST.resize(0);
  const float eps = FLT_EPSILON;
  const float epsd = UT::Inverse(BA_VARIANCE_MAX_DEPTH, BA_WEIGHT_FEATURE, eps);
  const float epsdST = UT::Inverse(BA_VARIANCE_MAX_DEPTH_SLIDING_TRACK, BA_WEIGHT_FEATURE, eps);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION)) {
      continue;
    }
    const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    const int iX = static_cast<int>(iX2d.size()), iXST = static_cast<int>(iXST2dST.size());
    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    iKF2X[iKF] = iX;          iKF2XST[iKF] = iXST;
    iX2d.resize(iX + Nx, -1); iXST2dST.resize(iXST + static_cast<int>(KF.m_STs.size()), -1);
    int *ix2d = iX2d.data() + iX, *ixST2dST = iXST2dST.data() + iXST;
    for (int ix = 0; ix < Nx; ++ix) {
      if (uds[ix] & LBA_FLAG_TRACK_UPDATE_INFORMATION) {
        if (KF.m_Axs[ix].m_Sadd.m_a > epsd) {
          ix2d[ix] = Nd++;
        } else if (!(uds[ix] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
          ix2d[ix] = -2;
        }
      }
      const int iST1 = KF.m_ix2ST[ix], iST2 = KF.m_ix2ST[ix + 1];
      if (iST2 == iST1) {
        continue;
      }
      ubyte update = LBA_FLAG_MARGINALIZATION_DEFAULT, nonZero = LBA_FLAG_MARGINALIZATION_DEFAULT;
      for (int iST = iST1; iST < iST2; ++iST) {
        if (KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION) {
          if (KF.m_AxsST[iST].m_Sadd.m_a > epsdST) {
            ixST2dST[iST] = NdST++;
            update = LBA_FLAG_MARGINALIZATION_UPDATE;
            nonZero = LBA_FLAG_MARGINALIZATION_NON_ZERO;
          } else if (!(KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
            ixST2dST[iST] = -2;
            update = LBA_FLAG_MARGINALIZATION_UPDATE;
          }
        } else {
          if (!(KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
            nonZero = LBA_FLAG_MARGINALIZATION_NON_ZERO;
          }
        }
      }
      KF.m_ms[ix] |= update;
      if (nonZero) {
        KF.m_ms[ix] |= LBA_FLAG_MARGINALIZATION_NON_ZERO;
      } else {
        KF.m_ms[ix] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
      }
    }
  }
  const int N = Nd + NdST, NC = SIMD_FLOAT_CEIL(N);
  m_work.Resize(NC + NC + Nd * sizeof(xp128f) / sizeof(float));
  float *mdds = m_work.Data(), *mddsST = mdds + Nd;
  float *nds = m_work.Data() + NC, *ndsST = nds + Nd;
  xp128f *_mdds = (xp128f *) (m_work.Data() + NC + NC);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    const KeyFrame &KF = m_KFs[iKF];
    const int *ix2d = iX2d.data() + iX;
    const int Nx = int(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      const int id = ix2d[ix];
      if (id < 0) {
        continue;
      }
      const FTR::Factor::DD &Sadd = KF.m_Axs[ix].m_Sadd;
      mdds[id] = Sadd.m_a;
      nds[id] = Sadd.m_b;
    }
    const int *ixST2dST = iXST2dST.data() + iKF2XST[iKF];
    const int NST = int(KF.m_STs.size());
    for (int iST = 0; iST < NST; ++iST) {
      const int idST = ixST2dST[iST];
      if (idST < 0) {
        continue;
      }
      const FTR::Factor::DD &SaddST = KF.m_AxsST[iST].m_Sadd;
      mddsST[idST] = SaddST.m_a;
      ndsST[idST] = SaddST.m_b;
    }
  }
  //SIMD::Add(N, UT::Inverse(BA_VARIANCE_REGULARIZATION_DEPTH, BA_WEIGHT_FEATURE), mdds);
  SIMD::Inverse(N, mdds);
  SIMD::Multiply(N, mdds, nds);
  for (int id = 0; id < Nd; ++id) {
    _mdds[id].vdup_all_lane(mdds[id]);
  }

  m_MxsTmp.Resize(NdST);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    const int iX = iKF2X[iKF];
    if (iX == -1) {
      continue;
    }
    const int iXST = iKF2XST[iKF];
    const int *ix2d = iX2d.data() + iX, *ixST2dST = iXST2dST.data() + iXST;
    const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      const int id = ix2d[ix];
      FTR::Factor::FixSource::Source::M &Mx = KF.m_Mxs[ix];
      if (id >= 0) {
        Mx.m_mdd.Set(mdds[id], nds[id]);
      }
      if (KF.m_Nsts[ix] == 0 || !(KF.m_ms[ix] & LBA_FLAG_MARGINALIZATION_UPDATE)) {
        continue;
      }
      const ubyte ud = uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH;
      const int iST1 = KF.m_ix2ST[ix], iST2 = KF.m_ix2ST[ix + 1];
      for (int iST = iST1; iST < iST2; ++iST) {
        const int idST = ixST2dST[iST];
        if (idST == -1) {
          continue;
        }
        FTR::Factor::FixSource::Source::M &MxST = KF.m_MxsST[iST];
        if (idST >= 0) {
          FTR::Factor::DD &dmddST = m_MxsTmp[idST].m_mdd;
          if (ud || (KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
            MxST.m_mdd.Set(mddsST[idST], ndsST[idST]);
            dmddST = MxST.m_mdd;
          } else {
            dmddST = MxST.m_mdd;
            MxST.m_mdd.Set(mddsST[idST], ndsST[idST]);
            FTR::Factor::DD::amb(MxST.m_mdd, dmddST, dmddST);
          }
        } else {
          if (!ud) {
            MxST.m_mdd.MakeMinus();
          }
        }
      }
    }
  }

  LA::ProductVector6f SmdczST;
  Camera::Factor::Unitary::CC dMczz;
  Camera::Factor::Binary::CC dMczm;
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    LocalFrame &LF = m_LFs[iLF];
    Camera::Factor::Unitary::CC &SMczz = m_SMcusLF[iLF];
    const bool ucz = (m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0;
    const int NZ = int(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      const int iKF = Z.m_iKF, iX = iKF2X[iKF];
      if (iX == -1) {
        continue;
      }
      const int *ix2d = iX2d.data() + iKF2X[iKF];
      const int *ixST2dST = iXST2dST.data() + iKF2XST[iKF];
      const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
      //const bool ucx = (m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0, ucr = ucx || ucz;
      const KeyFrame &KF = m_KFs[iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = LF.m_zs[iz].m_ix, id = ix2d[ix];
        if (id >= 0) {
          FTR::Factor::FixSource::Marginalize(_mdds[id], LF.m_Azs1[iz], &LF.m_Mzs1[iz]);
        }
        if (!(KF.m_ms[ix] & LBA_FLAG_MARGINALIZATION_UPDATE)) {
          continue;
        }
        FTR::Factor::DD &SmddST = LF.m_SmddsST[iz];
        const ubyte nonZero1 = (LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_NON_ZERO);
        if (KF.m_ms[ix] & LBA_FLAG_MARGINALIZATION_NON_ZERO) {
          const LocalFrame::SlidingTrack &ST = LF.m_STs[iz];
          const int iST0 = KF.m_ix2ST[ix], iST1 = iST0 + ST.m_ist1, iST2 = iST0 + ST.m_ist2;
          ubyte update = LBA_FLAG_MARGINALIZATION_DEFAULT;
          ubyte nonZero2 = LBA_FLAG_MARGINALIZATION_DEFAULT;
          for (int iST = iST1; iST < iST2; ++iST) {
            const int idST = ixST2dST[iST];
            if (idST >= 0) {
              if (ucz) {
                SmddST += KF.m_MxsST[iST].m_mdd;
              } else {
                SmddST += m_MxsTmp[idST].m_mdd;
              }
              update = LBA_FLAG_MARGINALIZATION_UPDATE;
              nonZero2 = LBA_FLAG_MARGINALIZATION_NON_ZERO;
            } else if (idST == -1) {
              if (!(KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
                if (ucz) {
                  SmddST += KF.m_MxsST[iST].m_mdd;
                }
                nonZero2 = LBA_FLAG_MARGINALIZATION_NON_ZERO;
              }
            } else {
              if (!ucz && !(uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH)) {
                SmddST += KF.m_MxsST[iST].m_mdd;
              }
              update = LBA_FLAG_MARGINALIZATION_UPDATE;
            }
          }
          LF.m_ms[iz] |= update;
          if (LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_UPDATE) {
            if (nonZero2) {
              if (!nonZero1) {
                LF.m_ms[iz] |= LBA_FLAG_MARGINALIZATION_NON_ZERO;
              }
              if (LF.m_Nsts[iz] == 1) {
                SmddST = KF.m_MxsST[iST1].m_mdd;
              }
            } else {
              if (nonZero1) {
                LF.m_ms[iz] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
                SmddST.MakeZero();
              }
            }
          }
        } else if (nonZero1) {
          LF.m_ms[iz] |= LBA_FLAG_MARGINALIZATION_UPDATE;
          LF.m_ms[iz] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
          SmddST.MakeZero();
        }
        if (!(LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_UPDATE)) {
          continue;
        }
        FTR::Factor::FixSource::M2 &M = LF.m_Mzs2[iz];
        if (LF.m_ms[iz] & LBA_FLAG_MARGINALIZATION_NON_ZERO) {
          if (!ucz && nonZero1) {
            dMczz = M.m_Mczz;
            FTR::Factor::FixSource::Marginalize(SmddST, LF.m_AzsST[iz].m_adcA, &SmdczST, &M);
            Camera::Factor::Unitary::CC::AmB(M.m_Mczz, dMczz, dMczz);
            SMczz += dMczz;
          } else {
            FTR::Factor::FixSource::Marginalize(SmddST, LF.m_AzsST[iz].m_adcA, &SmdczST, &M);
            SMczz += M.m_Mczz;
          }
        } else {
          if (!ucz && nonZero1) {
            M.m_Mczz.GetMinus(dMczz);
            SMczz += dMczz;
          }
        }

      }
    }

    const int NI = int(LF.m_Zm.m_Is.size());
    for (int iI = 0; iI < NI; ++iI) {
      const MeasurementMatchLF::Index &I = LF.m_Zm.m_Is[iI];
      const int iXST = iKF2XST[I.m_iKF];
      if (iXST == -1) {
        continue;
      }
      const int *ixST2dST = iXST2dST.data() + iXST;
      const ubyte *uds = m_uds.data() + m_iKF2d[I.m_iKF];
      Camera::Factor::Binary::CC &SMczm = LF.m_Zm.m_SMczms[I.m_ik];

      const int _iLF = LF.m_iLFsMatch[I.m_ik];
      const LocalFrame &_LF = m_LFs[_iLF];
      const bool _ucz = (m_ucsLF[_iLF] & LBA_FLAG_FRAME_UPDATE_CAMERA) != 0, uczm = ucz || _ucz;
      const KeyFrame &KF = m_KFs[I.m_iKF];
      const int i1 = LF.m_Zm.m_iI2zm[iI], i2 = LF.m_Zm.m_iI2zm[iI + 1];
      for (int i = i1; i < i2; ++i) {
        const FTR::Measurement::Match &izm = LF.m_Zm.m_izms[i];
        if (!(LF.m_ms[izm.m_iz1] & LBA_FLAG_MARGINALIZATION_UPDATE)) {
          continue;
        }
        float &SmddST = LF.m_Zm.m_SmddsST[i];
        const ubyte nonZero1 = (LF.m_Zm.m_ms[i] & LBA_FLAG_MARGINALIZATION_NON_ZERO);
        if (LF.m_ms[izm.m_iz1] & LBA_FLAG_MARGINALIZATION_NON_ZERO) {
          const int ix = LF.m_zs[izm.m_iz1].m_ix, iST0 = KF.m_ix2ST[ix];
          const int iST1 = iST0 + _LF.m_STs[izm.m_iz2].m_ist1;
          const int iST2 = iST0 + LF.m_STs[izm.m_iz1].m_ist2;
          ubyte update = LBA_FLAG_MARGINALIZATION_DEFAULT;
          ubyte nonZero2 = LBA_FLAG_MARGINALIZATION_DEFAULT;
          for (int iST = iST1; iST < iST2; ++iST) {
            const int idST = ixST2dST[iST];
            if (idST >= 0) {
              if (uczm) {
                SmddST += KF.m_MxsST[iST].m_mdd.m_a;
              } else {
                SmddST += m_MxsTmp[idST].m_mdd.m_a;
              }
              update = LBA_FLAG_MARGINALIZATION_UPDATE;
              nonZero2 = LBA_FLAG_MARGINALIZATION_NON_ZERO;
            } else if (idST == -1) {
              if (!(KF.m_usST[iST] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
                if (uczm) {
                  SmddST += KF.m_MxsST[iST].m_mdd.m_a;
                }
                nonZero2 = LBA_FLAG_MARGINALIZATION_NON_ZERO;
              }
            } else {
              if (!uczm && !(uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH)) {

                SmddST += KF.m_MxsST[iST].m_mdd.m_a;
              }
              update = LBA_FLAG_MARGINALIZATION_UPDATE;
            }
          }
          LF.m_Zm.m_ms[i] |= update;
          if (LF.m_Zm.m_ms[i] & LBA_FLAG_MARGINALIZATION_UPDATE) {
            if (nonZero2) {
              if (!nonZero1) {
                LF.m_Zm.m_ms[i] |= LBA_FLAG_MARGINALIZATION_NON_ZERO;
              }
              if (iST2 - iST1 == 1) {
                SmddST = KF.m_MxsST[iST1].m_mdd.m_a;
              }
            } else {
              if (nonZero1) {
                LF.m_Zm.m_ms[i] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
                SmddST = 0.0f;
              }
            }
          }
        } else if (nonZero1) {
          LF.m_Zm.m_ms[i] |= LBA_FLAG_MARGINALIZATION_UPDATE;
          LF.m_Zm.m_ms[i] &= ~LBA_FLAG_MARGINALIZATION_NON_ZERO;
          SmddST = 0.0f;
        }
        if (!(LF.m_Zm.m_ms[i] & LBA_FLAG_MARGINALIZATION_UPDATE)) {
          continue;
        }
        Camera::Factor::Binary::CC &Mczm = LF.m_Zm.m_Mczms[i];
        if (LF.m_Zm.m_ms[i] & LBA_FLAG_MARGINALIZATION_NON_ZERO) {
          if (!uczm && nonZero1) {
            dMczm = Mczm;
            FTR::Factor::FixSource::Marginalize(SmddST, LF.m_AzsST[izm.m_iz1],
                                                _LF.m_AzsST[izm.m_iz2], &SmdczST, &Mczm);
            Camera::Factor::Binary::CC::AmB(Mczm, dMczm, dMczm);
            SMczm += dMczm;
          } else {
            FTR::Factor::FixSource::Marginalize(SmddST, LF.m_AzsST[izm.m_iz1],
                                                _LF.m_AzsST[izm.m_iz2], &SmdczST, &Mczm);
            SMczm += Mczm;
          }

        } else {
          if (!uczm && nonZero1) {
            Mczm.GetMinus(dMczm);
            SMczm += dMczm;
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
    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (ix2d[ix] >= 0) {
        uds[ix] &= ~LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      } else if (ix2d[ix] == -2) {
        uds[ix] |= LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      }
      KF.m_ms[ix] &= ~LBA_FLAG_MARGINALIZATION_UPDATE;
    }
    const int *ixST2dST = iXST2dST.data() + iKF2XST[iKF];
    const int NST = static_cast<int>(KF.m_STs.size());
    for (int iST = 0; iST < NST; ++iST) {
      if (ixST2dST[iST] >= 0) {
        KF.m_usST[iST] &= ~LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      } else if (ixST2dST[iST] == -2) {
        KF.m_usST[iST] |= LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO;
      }
    }
  }
  for (int iLF = 0; iLF < nLFs; ++iLF) {
    LocalFrame &LF = m_LFs[iLF];
    const int NZ = int(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      if (iKF2X[Z.m_iKF] == -1) {
        continue;
      }
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        LF.m_ms[iz] &= ~LBA_FLAG_MARGINALIZATION_UPDATE;
      }
    }
    const int Nzm = int(LF.m_Zm.m_ms.size());
    for (int i = 0; i < Nzm; ++i) {
      LF.m_Zm.m_ms[i] &= ~LBA_FLAG_MARGINALIZATION_UPDATE;
    }
  }

}


bool LocalBundleAdjustor::SolveSchurComplement() {
  if (LBA_PROPAGATE_CAMERA >= 2 && m_iIter == 0 && SolveSchurComplementLast()) {
    return true;
  }
  bool scc = SolveSchurComplementPCG();

  if (LBA_EMBEDDED_MOTION_ITERATION) {
    EmbeddedMotionIteration();
  }
  if (!scc) {
    return false;
  }
  return true;
}

bool LocalBundleAdjustor::SolveSchurComplementPCG() {
  Camera::Factor::Unitary::CC Acc;
  const int pc = 6, pm = 9;
  const int pcm = pc + pm;
  const int Nc = int(m_LFs.size()), Ncp = Nc * pc, Nmp = Nc * pm, N = Ncp + Nmp;
  m_Acus.Resize(Nc);
  m_Amus.Resize(Nc);
  m_bs.Resize(N);
  float *bc = m_bs.Data(), *bm = bc + Ncp;
  const int Nb = CountSchurComplementsOffDiagonal();
  m_Acbs.Resize(Nb);
  m_ic2b.resize(Nc);
  const float ar = UT::Inverse(BA_VARIANCE_REGULARIZATION_ROTATION, BA_WEIGHT_FEATURE);
  const float ap = UT::Inverse(BA_VARIANCE_REGULARIZATION_POSITION, BA_WEIGHT_FEATURE);
  const float av = UT::Inverse(BA_VARIANCE_REGULARIZATION_VELOCITY, BA_WEIGHT_FEATURE);
  const float aba = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_ACCELERATION, BA_WEIGHT_FEATURE);
  const float abw = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE);
  for (int ic = 0, ib = 0; ic < Nc; ++ic, bc += pc, bm += pm) {
    m_ic2b[ic] = ib;
    const int iLF = m_ic2LF[ic];
    const LocalFrame &LF = m_LFs[iLF];
    Camera::Factor::Unitary::CC::AmB(m_SAcusLF[iLF], m_SMcusLF[iLF], Acc);
    Acc.m_A.IncreaseDiagonal(ap, ar);
    Acc.m_A.GetAlignedMatrix6x6f(m_Acus[ic]);
    Acc.m_b.Get(bc);
    const Camera::Factor::Unitary::MM &Amm = m_SAcmsLF[iLF].m_Au.m_Amm;
    m_Amus[ic].Set(Amm.m_A);
    m_Amus[ic].IncreaseDiagonal(av, aba, abw);
    Amm.m_b.Get(bm);
    const int _ic = ic + 1;
    if (_ic == Nc) {
      continue;
    }
    LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + ib;
    const int _iLF = m_ic2LF[_ic];
    const Camera::Factor::Binary &SAcmb = m_SAcmsLF[_iLF].m_Ab;
    LA::AlignedMatrix6x6f::AmB(SAcmb.m_Acc, LF.m_Zm.m_SMczms[0], Acbs[0]);
    const int Nk = static_cast<int>(LF.m_iLFsMatch.size());

    for (int ik = 1; ik < Nk; ++ik) {
      LF.m_Zm.m_SMczms[ik].GetMinus(Acbs[ik]);
    }
    ib += Nk;
  }
  m_AcbTs.Resize(Nb);
  for (int ib = 0; ib < Nb; ++ib) {
    m_Acbs[ib].GetTranspose(m_AcbTs[ib]);
  }
  PrepareConditioner();

  m_xs.Resize(N);
  m_rs.Resize(N);
  m_ps.Resize(N);
  m_zs.Resize(N);
  m_drs.Resize(N);
  m_dxs.Resize(N);

  bool scc = true;
  float Se2, Se2Pre, Se2Min, e2Max, e2MaxMin, alpha, beta;
  m_rs = m_bs;
#ifdef CFG_INCREMENTAL_PCG
  LA::Vector6f *xcs = (LA::Vector6f *) m_xs.Data();
  LA::Vector9f *xms = (LA::Vector9f *) (xcs + Nc);
  for (int ic = 0; ic < Nc; ++ic) {
    const int iLF = m_ic2LF[ic];
    xcs[ic] = m_xcsLF[iLF];
    xms[ic] = m_xmsLF[iLF];
  }
  m_xs.MakeMinus();
  m_xsGN = m_xs;
  ApplyA(m_xs, &m_drs);
  m_rs -= m_drs;
#else
  m_xsGN.Resize(N);
  m_xsGN.MakeZero();
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
    const float e2MaxConv[2] = {ME::ChiSquareDistance<3>(BA_PCG_MIN_CONVERGE_PROBABILITY,
                                                         BA_WEIGHT_FEATURE),
                                ME::ChiSquareDistance<3>(BA_PCG_MAX_CONVERGE_PROBABILITY,
                                                         BA_WEIGHT_FEATURE)};
    Se2Min = Se2;
    e2MaxMin = e2Max;

    m_drs *= alpha;
    m_rs -= m_drs;
#ifdef CFG_INCREMENTAL_PCG
    m_ps.GetScaled(alpha, m_dxs);
    m_xs += m_dxs;
#else
    m_ps.GetScaled(alpha, m_xs);
#endif
#ifdef CFG_INCREMENTAL_PCG
    ApplyM(m_bs, &m_drs);
    Se2Pre = m_bs.Dot(m_drs);
    const float Se2ConvMin = Se2Pre * BA_PCG_MIN_CONVERGE_RESIDUAL_RATIO;
    const float Se2ConvMax = Se2Pre * BA_PCG_MAX_CONVERGE_RESIDUAL_RATIO;
#else
    const float Se2ConvMin = Se2 * BA_PCG_MIN_CONVERGE_RESIDUAL_RATIO;
    const float Se2ConvMax = Se2 * BA_PCG_MAX_CONVERGE_RESIDUAL_RATIO;
#endif
    int cnt = 0;
    const int nIters = std::min(N, BA_PCG_MAX_ITERATIONS);
    for (m_iIterPCG = 0; m_iIterPCG < nIters; ++m_iIterPCG) {
      ApplyM(m_rs, &m_zs);
      Se2Pre = Se2;
      ConvertCameraMotionResiduals(m_rs, m_zs, &Se2, &e2Max);

      if (Se2 < Se2Min) {
        Se2Min = Se2;
        e2MaxMin = e2Max;
        m_xsGN = m_xs;
        //cnt = 0;
      } else {
        ++cnt;
      }

      const int i = (Se2Min <= Se2ConvMin && m_iIterPCG >= BA_PCG_MIN_ITERATIONS) ? 0 : 1;
      if (Se2 == 0.0f || e2MaxMin < e2MaxConv[i]) {
        scc = true;
        break;
      } else if (Se2Min > Se2ConvMax) {
        scc = false;
        break;
      }
      beta = Se2 / Se2Pre;

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


bool LocalBundleAdjustor::SolveSchurComplementLast() {
  const int Nc = static_cast<int>(m_LFs.size());
  if (Nc < 3) {
    return false;
  }
  const int pc = 6, pm = 9;
  const int pcm = pc + pm/*, pmcm = pm + pcm*/;
  const float eps = FLT_EPSILON;
  const float epsr = UT::Inverse(BA_VARIANCE_MAX_ROTATION, BA_WEIGHT_FEATURE, eps);
  const float epsp = UT::Inverse(BA_VARIANCE_MAX_POSITION, BA_WEIGHT_FEATURE, eps);
  const float epsv = UT::Inverse(BA_VARIANCE_MAX_VELOCITY, BA_WEIGHT_FEATURE, eps);
  const float epsba = UT::Inverse(BA_VARIANCE_MAX_BIAS_ACCELERATION, BA_WEIGHT_FEATURE, eps);
  const float epsbw = UT::Inverse(BA_VARIANCE_MAX_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE, eps);

  const float _eps[pcm] = {epsp, epsp, epsp, epsr, epsr, epsr,
                           epsv, epsv, epsv, epsba, epsba, epsba, epsbw, epsbw, epsbw};
  const float ar = UT::Inverse(BA_VARIANCE_REGULARIZATION_ROTATION, BA_WEIGHT_FEATURE);
  const float ap = UT::Inverse(BA_VARIANCE_REGULARIZATION_POSITION, BA_WEIGHT_FEATURE);
  const float av = UT::Inverse(BA_VARIANCE_REGULARIZATION_VELOCITY, BA_WEIGHT_FEATURE);
  const float aba = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_ACCELERATION, BA_WEIGHT_FEATURE);
  const float abw = UT::Inverse(BA_VARIANCE_REGULARIZATION_BIAS_GYROSCOPE, BA_WEIGHT_FEATURE);
  
  LA::AlignedMatrixXf A;
  LA::AlignedVectorXf b;

  m_work.Resize(A.BindSize(pcm, pcm) + b.BindSize(pcm));
  A.Bind(m_work.Data(), pcm, pcm);
  b.Bind(A.BindNext(), pcm);

  Camera::Factor::Unitary::CC Acc;
  LA::AlignedMatrix9x9f Amm;

  const int ic2 = Nc - 1, iLF2 = m_ic2LF[ic2];
  Camera::Factor::Unitary::CC::AmB(m_SAcusLF[iLF2], m_SMcusLF[iLF2], Acc);
  Acc.m_A.IncreaseDiagonal(ap, ar);
  A.SetBlockDiagonal(0, Acc.m_A);
  const Camera::Factor &Acm2 = m_SAcmsLF[iLF2];
  A.SetBlock(0, pc, Acm2.m_Au.m_Acm);
  b.SetBlock(0, Acc.m_b);
  Amm.Set(Acm2.m_Au.m_Amm.m_A);
  Amm.IncreaseDiagonal(av, aba, abw);
  A.SetBlock(pc, pc, Amm);
  b.SetBlock(pc, Acm2.m_Au.m_Amm.m_b);

  A.SetLowerFromUpper();
  if (!A.SolveLDL(b, _eps)) {
    return false;
  }
  b.MakeMinus();

  LA::AlignedVector3f x;
  const int N = Nc * pcm;
  m_xsGN.Resize(N);     m_xsGN.MakeZero();
  m_xp2s.Resize(Nc);    m_xp2s.MakeZero();
  m_xr2s.Resize(Nc);    m_xr2s.MakeZero();
  m_xv2s.Resize(Nc);    m_xv2s.MakeZero();
  m_xba2s.Resize(Nc);   m_xba2s.MakeZero();
  m_xbw2s.Resize(Nc);   m_xbw2s.MakeZero();
  float *xc2 = m_xsGN.Data() + ic2 * pc;
  float *xm2 = m_xsGN.End() - pm;
  b.GetBlock(0, x);   x.Get(xc2);     m_xp2s[ic2] = x.SquaredLength();
  b.GetBlock(3, x);   x.Get(xc2 + 3); m_xr2s[ic2] = x.SquaredLength();
  b.GetBlock(6, x);   x.Get(xm2);     m_xv2s[ic2] = x.SquaredLength();
  b.GetBlock(9, x);   x.Get(xm2 + 3); m_xba2s[ic2] = x.SquaredLength();
  b.GetBlock(12, x);  x.Get(xm2 + 6); m_xbw2s[ic2] = x.SquaredLength();
  return true;
}

void LocalBundleAdjustor::PrepareConditioner() {
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
  const int Nb = std::min(LBA_PCG_CONDITIONER_BAND, LBA_MAX_SLIDING_TRACK_LENGTH);
  const int Nc = static_cast<int>(m_LFs.size());
  if (Nb <= 1) {
    m_Mcs.Resize(Nc);
    m_Mms.Resize(Nc);
    for (int ic = 0; ic < Nc; ++ic) {
      const int iLF = m_ic2LF[ic];
      m_Mcs[ic].Set(m_SAcusLF[iLF], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsc);
      //m_Mcs[ic].Set(m_Acus[ic], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsc);
      m_Mms[ic].Set(m_Amus[ic], BA_PCG_CONDITIONER_MAX, BA_PCG_CONDITIONER_EPSILON, epsm);
    }
    return;
  }
  m_Mcc.Resize(Nc, Nb);     m_MccT.Resize(Nc, Nb);
  m_Mcm.Resize(Nc, 2);      m_McmT.Resize(Nc, 2);
  m_Mmc.Resize(Nc, Nb - 1); m_MmcT.Resize(Nc, Nb - 1);
  m_Mmm.Resize(Nc, 2);      m_MmmT.Resize(Nc, 2);
  for (int ic = 0; ic < Nc; ++ic) {
    LA::AlignedMatrix6x6f *Accs = m_Mcc[ic];
    Accs[0] = m_Acus[ic];
    const LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + m_ic2b[ic] - 1;
    const int Nbc = ic + Nb > Nc ? Nc - ic : Nb;
    for (int ib = 1; ib < Nbc; ++ib) {
      Accs[ib] = Acbs[ib];
    }
    m_Mcm[ic][0] = m_SAcmsLF[m_ic2LF[ic]].m_Au.m_Acm;
    m_Mmm[ic][0] = m_Amus[ic];
    const int _ic = ic + 1;
    if (_ic == Nc) {
      continue;
    }
    const Camera::Factor::Binary &Ab = m_SAcmsLF[m_ic2LF[_ic]].m_Ab;
    m_Mcm[ic][1] = Ab.m_Acm;
    m_Mmm[ic][1] = Ab.m_Amm;
    LA::AlignedMatrix9x6f *Amcs = m_Mmc[ic];
    Amcs[0] = Ab.m_Amc;
    const int Nbm = Nbc - 1;
    for (int ib = 1; ib < Nbm; ++ib) {
      Amcs[ib].MakeZero();
    }
  }

  AlignedVector<LA::AlignedMatrix6x6f> AccsT;
  AlignedVector<LA::AlignedMatrix9x6f> AcmsT;
  AlignedVector<LA::AlignedMatrix6x9f> AmcsT;
  AlignedVector<LA::AlignedMatrix9x9f> AmmsT;
  m_work.Resize((AccsT.BindSize(Nb) + AcmsT.BindSize(2) +
                 AmcsT.BindSize(Nb) + AmmsT.BindSize(2)) / sizeof(float));
  AccsT.Bind(m_work.Data(), Nb);
  AcmsT.Bind(AccsT.BindNext(), 2);
  AmcsT.Bind(AcmsT.BindNext(), Nb);
  AmmsT.Bind(AmcsT.BindNext(), Nb);
  for (int ic = 0; ic < Nc; ++ic) {
    LA::AlignedMatrix6x6f *Mccs = m_Mcc[ic], *MccsT = m_MccT[ic];
    LA::AlignedMatrix6x9f *Mcms = m_Mcm[ic], *MmcsT = m_MmcT[ic];
    LA::AlignedMatrix9x6f *McmsT = m_McmT[ic], *Mmcs = m_Mmc[ic];
    LA::AlignedMatrix9x9f *Mmms = m_Mmm[ic], *MmmsT = m_MmmT[ic];
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = ic + 1 == Nc ? 1 : 2;
    const int Nbmc = Nbcc - 1;
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
        Mcms[ib].GetTranspose(AcmsT[ib]);
        LA::AlignedMatrix9x6f::ABT(Mcc, AcmsT[ib], Mcms[ib]);
        Mcms[ib].GetTranspose(McmsT[ib]);
      }
      for (int ib = 1; ib < Nbcc; ++ib) {
        const LA::AlignedMatrix6x6f &MccT = MccsT[ib];
        const int _ic = ic + ib;
        LA::AlignedMatrix6x6f *_Mccs = m_Mcc[_ic] - ib;
        LA::AlignedMatrix6x6f::AddABTToUpper(MccT, AccsT[ib], _Mccs[ib]);

        for (int jb = ib + 1; jb < Nbcc; ++jb) {
          LA::AlignedMatrix6x6f::AddABTTo(MccT, AccsT[jb], _Mccs[jb]);
        }
        if (ib == 1) {
          LA::AlignedMatrix9x6f::AddABTTo(MccT, AcmsT[ib], m_Mcm[_ic][0]);
        }
      }
      for (int ib = 0; ib < Nbcm; ++ib) {
        const LA::AlignedMatrix9x6f &McmT = McmsT[ib];
        const int _ic = ic + ib;
        const LA::AlignedMatrix6x6f *_AccsT = AccsT.Data() + 1;
        LA::AlignedMatrix9x6f *_Mmcs = m_Mmc[_ic] - ib;
        for (int jb = ib; jb < Nbmc; ++jb) {
          LA::AlignedMatrix9x6f::AddABTTo(McmT, _AccsT[jb], _Mmcs[jb]);
        }
        LA::AlignedMatrix9x9f *_Mmms = m_Mmm[_ic];
        LA::AlignedMatrix9x9f::AddABTToUpper(McmT, AcmsT[ib], _Mmms[0]);

        if (ib == 0 && Nbcm == 2) {
          LA::AlignedMatrix9x9f::AddABTTo(McmT, AcmsT[1], _Mmms[1]);
        }
      }
    } else {
      for (int ib = 0; ib < Nbcc; ++ib) {
        Mccs[ib].MakeZero();
        MccsT[ib].MakeZero();
      }
      for (int ib = 0; ib < Nbcm; ++ib) {
        Mcms[ib].MakeZero();
        McmsT[ib].MakeZero();
      }
    }
    LA::AlignedMatrix9x9f &Mmm = Mmms[0];
    if (Mmm.InverseLDL(epsm)) {
      MmmsT[0] = Mmm;
      Mmm.MakeMinus();
      for (int ib = 0; ib < Nbmc; ++ib) {
        Mmcs[ib].GetTranspose(AmcsT[ib]);
        LA::AlignedMatrix9x9f::ABT(Mmm, AmcsT[ib], Mmcs[ib]);
        Mmcs[ib].GetTranspose(MmcsT[ib]);
      }
      if (Nbmm == 2) {
        Mmms[1].GetTranspose(AmmsT[1]);
        LA::AlignedMatrix9x9f::ABT(Mmm, AmmsT[1], Mmms[1]);
        Mmms[1].GetTranspose(MmmsT[1]);
      }
      for (int ib = 0; ib < Nbmc; ++ib) {
        const LA::AlignedMatrix6x9f &MmcT = MmcsT[ib];
        const int _ic = ic + ib + 1;
        LA::AlignedMatrix6x6f *_Mccs = m_Mcc[_ic] - ib;
        LA::AlignedMatrix6x9f::AddABTToUpper(MmcT, AmcsT[ib], _Mccs[ib]);

        for (int jb = ib + 1; jb < Nbmc; ++jb) {
          LA::AlignedMatrix6x9f::AddABTTo(MmcT, AmcsT[jb], _Mccs[jb]);
        }
        if (ib == 0 && Nbmm == 2) {
          LA::AlignedMatrix9x9f::AddABTTo(MmcT, AmmsT[1], m_Mcm[_ic][0]);
        }
      }
      if (Nbmm == 2) {
        const LA::AlignedMatrix9x9f &MmmT = MmmsT[1];
        const int _ic = ic + 1;
        LA::AlignedMatrix9x6f *_Mmcs = m_Mmc[_ic] - 1;
        for (int jb = 1; jb < Nbmc; ++jb) {
          LA::AlignedMatrix9x9f::AddABTTo(MmmT, AmcsT[jb], _Mmcs[jb]);
        }
        LA::AlignedMatrix9x9f::AddABTToUpper(MmmT, AmmsT[1], m_Mmm[_ic][0]);

      }
    } else {
      for (int ib = 0; ib < Nbmc; ++ib) {
        Mmcs[ib].MakeZero();
        MmcsT[ib].MakeZero();
      }
      for (int ib = 0; ib < Nbmm; ++ib) {
        Mmms[ib].MakeZero();
        MmmsT[ib].MakeZero();
      }
    }
  }
}

void LocalBundleAdjustor::ApplyM(const LA::AlignedVectorXf &xs, LA::AlignedVectorXf *Mxs) {
  const int Nb = std::min(LBA_PCG_CONDITIONER_BAND, LBA_MAX_SLIDING_TRACK_LENGTH);
  const int Nc = int(m_LFs.size());
  if (Nb <= 1) {
#ifdef CFG_PCG_FULL_BLOCK
    LA::ProductVector6f xc;
    LA::AlignedVector9f xm;
#else
    LA::AlignedVector3f xp, xr, xv, xba, xbw;
#endif
    Mxs->Resize(xs.Size());
    const LA::Vector6f *xcs = (LA::Vector6f *) xs.Data();
    LA::Vector6f *Mxcs = (LA::Vector6f *) Mxs->Data();
    for (int ic = 0; ic < Nc; ++ic) {
#ifdef CFG_PCG_FULL_BLOCK
      xc.Set(xcs[ic]);
      m_Mcs[ic].Apply(xc, (PCG_TYPE *) &Mxcs[ic]);
#else
      xcs[ic].Get(xp, xr);
      m_Mcs[ic].Apply(xp, xr, (float *) Mxcs[ic]);
#endif
    }
    const LA::Vector9f *xms = (LA::Vector9f *) (xcs + Nc);
    LA::Vector9f *Mxms = (LA::Vector9f *) (Mxcs + Nc);
    for (int ic = 0; ic < Nc; ++ic) {
#ifdef CFG_PCG_FULL_BLOCK
      xm.Set(xms[ic]);
      m_Mms[ic].Apply(xm, (PCG_TYPE *) &Mxms[ic]);
#else
      xms[ic].Get(xv, xba, xbw);
      m_Mms[ic].Apply(xv, xba, xbw, (float *) Mxms[ic]);
#endif
    }
    return;
  }
  LA::ProductVector6f bc;
  LA::AlignedVector9f bm;
  Mxs->Set(xs);
  LA::Vector6f *bcs = (LA::Vector6f *) Mxs->Data();
  LA::Vector9f *bms = (LA::Vector9f *) (bcs + Nc);
  for (int ic = 0; ic < Nc; ++ic) {
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = ic + 1 == Nc ? 1 : 2;
    const int Nbmc = Nbcc - 1;
    const int Nbmm = Nbcm;
    
    bc.Set(bcs[ic]);
    const LA::AlignedMatrix6x6f *MccsT = m_MccT[ic];
    for (int ib = 1; ib < Nbcc; ++ib) {
      LA::AlignedMatrix6x6f::AddAbTo<float>(MccsT[ib], bc, bcs[ic + ib]);
    }
    const LA::AlignedMatrix9x6f *McmsT = m_McmT[ic];
    for (int ib = 0; ib < Nbcm; ++ib) {
      LA::AlignedMatrix9x6f::AddAbTo<float>(McmsT[ib], bc, bms[ic + ib]);
    }
    LA::AlignedMatrix6x6f::Ab<float>(MccsT[0], bc, bcs[ic]);

    bm.Set(bms[ic]);
    const LA::AlignedMatrix6x9f *MmcsT = m_MmcT[ic];
    for (int ib = 0; ib < Nbmc; ++ib) {
      LA::AlignedMatrix6x9f::AddAbTo<float>(MmcsT[ib], bm, bcs[ic + ib + 1]);
    }
    const LA::AlignedMatrix9x9f *MmmsT = m_MmmT[ic];
    if (Nbmm == 2) {
      LA::AlignedMatrix9x9f::AddAbTo<float>(MmmsT[1], bm, bms[ic + 1]);
    }
    LA::AlignedMatrix9x9f::Ab<float>(MmmsT[0], bm, bms[ic]);
  }
  m_bcs.Resize(Nc);
  m_bms.Resize(Nc);
  for (int ic = Nc - 1; ic >= 0; --ic) {
    const int Nbcc = ic + Nb > Nc ? Nc - ic : Nb;
    const int Nbcm = ic + 1 == Nc ? 1 : 2;
    const int Nbmc = Nbcc - 1;
    const int Nbmm = Nbcm;

    float *_bm = bms[ic];
    const LA::AlignedMatrix9x6f *Mmcs = m_Mmc[ic];
    for (int ib = 0; ib < Nbmc; ++ib) {
      LA::AlignedMatrix9x6f::AddAbTo(Mmcs[ib], m_bcs[ic + ib + 1], _bm);
    }
    if (Nbmm == 2) {
      LA::AlignedMatrix9x9f::AddAbTo(m_Mmm[ic][1], m_bms[ic + 1], _bm);
    }
    m_bms[ic].Set(_bm);

    float *_bc = bcs[ic];
    const LA::AlignedMatrix6x6f *Mccs = m_Mcc[ic];
    for (int ib = 1; ib < Nbcc; ++ib) {
      LA::AlignedMatrix6x6f::AddAbTo(Mccs[ib], m_bcs[ic + ib], _bc);
    }
    const LA::AlignedMatrix6x9f *Mcms = m_Mcm[ic];
    for (int ib = 0; ib < Nbcm; ++ib) {
      LA::AlignedMatrix6x9f::AddAbTo(Mcms[ib], m_bms[ic + ib], _bc);
    }
    m_bcs[ic].Set(_bc);
  }
}

void LocalBundleAdjustor::ApplyA(const LA::AlignedVectorXf &xs, LA::AlignedVectorXf *Axs) {
  //LA::AlignedVector6f v6;
  const LA::Vector6f *xcs = (LA::Vector6f *) xs.Data();
  ConvertCameraUpdates(xcs, &m_xcsP);
  Axs->Resize(xs.Size());
  LA::Vector6f *Axcs = (LA::Vector6f *) Axs->Data();
  const int Nc = int(m_LFs.size());
  if (m_Acus.Size() < Nc) {
    Axs->MakeZero();
    return;
  }
  //m_Axcs.Resize(Nc);
  for (int ic = 0; ic < Nc; ++ic) {
    LA::AlignedMatrix6x6f::Ab(m_Acus[ic], m_xcsP[ic], (float *) &Axcs[ic]);
  }
  for (int ic = 0; ic < Nc; ++ic) {
    const int ib = m_ic2b[ic];
    const LA::AlignedMatrix6x6f *Acbs = m_Acbs.Data() + ib;
    const LA::AlignedMatrix6x6f *AcbTs = m_AcbTs.Data() + ib;
    const LA::ProductVector6f &xc = m_xcsP[ic];
    float *Axc = Axcs[ic];
    const LocalFrame &LF = m_LFs[m_ic2LF[ic]];
    const int Nk = int(LF.m_iLFsMatch.size());
    for (int ik = 0; ik < Nk; ++ik) {
      const int _ic = ic + ik + 1;
      LA::AlignedMatrix6x6f::AddAbTo(Acbs[ik], m_xcsP[_ic], Axc);
      LA::AlignedMatrix6x6f::AddAbTo(AcbTs[ik], xc, (float *) &Axcs[_ic]);
    }
  }
  //ConvertCameraUpdates(m_Axcs, Axcs);
  ApplyAcm(m_xcsP.Data(), (LA::Vector9f *) (xcs + Nc), Axcs, (LA::Vector9f *) (Axcs + Nc), false,
           m_Amus.Size() == Nc ? m_Amus.Data() : NULL);
}

void LocalBundleAdjustor::ApplyAcm(const LA::ProductVector6f *xcs, const LA::Vector9f *xms,
                                   LA::Vector6f *Axcs, LA::Vector9f *Axms, const bool Acc,
                                   const LA::AlignedMatrix9x9f *Amus) {
  LA::AlignedMatrix6x6f A66;
  LA::AlignedMatrix6x9f A69;
  LA::AlignedMatrix9x6f A96;
  LA::AlignedMatrix9x9f A99;
  LA::AlignedVector9f v9[2];

  LA::AlignedMatrix3x3f A33;
  LA::AlignedMatrix3x9f A39;
  LA::AlignedVector3f v3;

  const int Nc = int(m_LFs.size());
  for (int ic = 0, r = 0; ic < Nc; ++ic, r = 1 - r) {
    const LA::ProductVector6f &xc = xcs[ic];
    LA::AlignedVector9f &xm = v9[r];
    xm.Set(xms[ic]);
    float *Axc = Axcs[ic], *Axm = Axms[ic];
    const Camera::Factor &SAcm = m_SAcmsLF[m_ic2LF[ic]];
    if (Amus) {
      LA::AlignedMatrix9x9f::Ab(Amus[ic], xm, Axm);
    } else {
      A99.Set(SAcm.m_Au.m_Amm.m_A);
      LA::AlignedMatrix9x9f::Ab(A99, xm, Axm);
    }
    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Au.m_Acm, xm, Axc);
    SAcm.m_Au.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, xc, Axm);
    if (ic == 0) {
      continue;
    }
    const int _ic = ic - 1;
    const LA::ProductVector6f &_xc = xcs[_ic];
    const LA::AlignedVector9f &_xm = v9[1 - r];
    float *_Axc = Axcs[_ic], *_Axm = Axms[_ic];
    if (Acc) {
      LA::AlignedMatrix6x6f::AddAbTo(SAcm.m_Ab.m_Acc, xc, _Axc);
      SAcm.m_Ab.m_Acc.GetTranspose(A66);
      LA::AlignedMatrix6x6f::AddAbTo(A66, _xc, Axc);
    }

    LA::AlignedMatrix6x9f::AddAbTo(SAcm.m_Ab.m_Acm, xm, _Axc);
    SAcm.m_Ab.m_Acm.GetTranspose(A96);
    LA::AlignedMatrix9x6f::AddAbTo(A96, _xc, Axm);
    LA::AlignedMatrix9x6f::AddAbTo(SAcm.m_Ab.m_Amc, xc, _Axm);
    SAcm.m_Ab.m_Amc.GetTranspose(A69);
    LA::AlignedMatrix6x9f::AddAbTo(A69, _xm, Axc);
    LA::AlignedMatrix9x9f::AddAbTo(SAcm.m_Ab.m_Amm, xm, _Axm);
    SAcm.m_Ab.m_Amm.GetTranspose(A99);
    LA::AlignedMatrix9x9f::AddAbTo(A99, _xm, Axm);
  }
}

LocalBundleAdjustor::Residual LocalBundleAdjustor::ComputeResidual(const LA::AlignedVectorXf &xs,
                                                                   const bool minus) {
  Residual R;
  ApplyA(xs, &m_rs);
  if (minus) {
    R.m_F = xs.Dot(m_rs) / 2 - xs.Dot(m_bs);
    m_rs -= m_bs;
  } else {
    R.m_F = xs.Dot(m_rs) / 2 + xs.Dot(m_bs);
    m_rs += m_bs;
  }
  R.m_r2 = m_rs.SquaredLength();
  return R;
}

void LocalBundleAdjustor::SolveBackSubstitution() {
  const int nKFs = static_cast<int>(m_KFs.size());
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION)) {
      continue;
    }
    const KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    if (Nx == 0) {
      continue;
    }
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & LBA_FLAG_TRACK_UPDATE_INFORMATION)) {
        continue;
      }
      uds[ix] |= LBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION;
      m_ucsKF[iKF] |= LBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION;
    }
  }
  const int Nc = static_cast<int>(m_LFs.size());
  for (int ic = 0; ic < Nc; ++ic) {
    const int iLF = m_ic2LF[ic];
    if (m_xr2s[ic] <= BA_BACK_SUBSTITUTE_ROTATION &&
        m_xp2s[ic] <= BA_BACK_SUBSTITUTE_POSITION) {
      m_ucsLF[iLF] &= ~LBA_FLAG_FRAME_UPDATE_DELTA;
      continue;
    }
    m_ucsLF[iLF] |= LBA_FLAG_FRAME_UPDATE_DELTA;
    const LocalFrame &LF = m_LFs[iLF];
    const int NZ = int(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      if (Z.m_iz1 == Z.m_iz2) {
        continue;
      }
      m_ucsKF[Z.m_iKF] |= LBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION;
      ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        uds[LF.m_zs[iz].m_ix] |= LBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION;
      }
    }
  }

  int iX = 0;
  m_iKF2X.assign(nKFs, -1);
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION)) {
      continue;
    }
    m_iKF2X[iKF] = iX;
    iX += static_cast<int>(m_KFs[iKF].m_xs.size());
  }
  m_xds.Resize(iX);

  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_BACK_SUBSTITUTION)) {
      continue;
    }
    const ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    float *xds = m_xds.Data() + m_iKF2X[iKF];
    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      if (!(uds[ix] & LBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION)) {
        continue;
      } else if (uds[ix] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO) {
        xds[ix] = 0.0f;
      } else {
        xds[ix] = KF.m_Mxs[ix].BackSubstitute();
      }
    }
  }
  LA::AlignedVector6f xc;
  const LA::Vector6f *xcs = (LA::Vector6f *) m_xsGN.Data();
  for (int ic = 0; ic < Nc; ++ic) {
    const int iLF = m_ic2LF[ic];
    if (!(m_ucsLF[iLF] & LBA_FLAG_FRAME_UPDATE_DELTA)) {
      continue;
    }
    xc.Set(xcs[ic]);
    const LocalFrame &LF = m_LFs[iLF];
    const int NZ = int(LF.m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = LF.m_Zs[iZ];
      const ubyte *uds = m_uds.data() + m_iKF2d[Z.m_iKF];
      float *xds = m_xds.Data() + m_iKF2X[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const int ix = LF.m_zs[iz].m_ix;
        if (!(uds[ix] & LBA_FLAG_TRACK_UPDATE_INFORMATION_ZERO)) {
          xds[ix] = LF.m_Mzs1[iz].BackSubstitute(xc) + xds[ix];
        }
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
      if (!(uds[ix] & LBA_FLAG_TRACK_UPDATE_BACK_SUBSTITUTION)) {
        continue;
      }
      xds[ix] = -xds[ix];
    }
  }
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(m_ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION)) {
      continue;
    }
    m_ucsKF[iKF] &= ~LBA_FLAG_FRAME_UPDATE_TRACK_INFORMATION;
    ubyte *uds = m_uds.data() + m_iKF2d[iKF];
    KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
    for (int ix = 0; ix < Nx; ++ix) {
      uds[ix] &= ~LBA_FLAG_TRACK_UPDATE_INFORMATION;
    }
    const int NST = static_cast<int>(KF.m_STs.size());
    for (int iST = 0; iST < NST; ++iST) {
      KF.m_usST[iST] &= ~LBA_FLAG_TRACK_UPDATE_INFORMATION;
    }
  }

  PushDepthUpdates(m_xds, &m_xsGN);
  m_x2GN = m_xsGN.SquaredLength();
}

bool LocalBundleAdjustor::EmbeddedMotionIteration() {
  const int pc = 6, pm = 9;
  const int Nc = m_CsLF.Size();
  const LA::Vector6f *xcs = (LA::Vector6f *) m_xsGN.Data();
  LA::Vector9f *xms = (LA::Vector9f *) (xcs + Nc);
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
  for (int ic1 = -1, ic2 = 0, imp1 = -pm, imp2 = 0, r = 0; ic2 < Nc;
       ic1 = ic2++, imp1 = imp2, imp2 += pm, r = 1 - r) {
    const Camera::Factor &_A = m_SAcmsLF[m_ic2LF[ic2]];
    Amm.Set(_A.m_Au.m_Amm.m_A);
    A.SetBlock(imp2, 0, Amm);
    float *_b = b.Data() + imp2;
    _A.m_Au.m_Amm.m_b.Get(_b);
    xc[r].Set(xcs[ic2]);
    _A.m_Au.m_Acm.GetTranspose(Amc);
    LA::AlignedMatrix9x6f::AddAbTo(Amc, xc[r], _b);
    if (ic2 == 0) {
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
  for (int ic = 0, imp = 0; ic < Nc; ++ic) {
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
      if (ic == Nc - 1) {
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

  for (int ic = 0, imp = 0; ic < Nc; ++ic, imp += pm) {
    memcpy(x[ic], b.Data() + imp, 36);
  }
  for (int ic = Nc - 1, imp = Nmr - 1, r = ic & 1; ic >= 0; --ic, r = 1 - r) {
    const int _ic = ic + 1;
    const int _Nmc = _ic == Nc ? pm : Nmc;
    float *xi = x[ic];
    if (_ic < Nc) {
      memcpy(xi + pm, x[_ic], 36);
    }
    for (int ip = pm - 1; ip >= 0; --ip, --imp) {
      xi[ip] -= SIMD::Dot(ip + 1, _Nmc, A[imp], xi);
    }
  }
  for (int ic = 0; ic < Nc; ++ic) {
    xms[ic].Set(x[ic]);
  }

  m_xsGN.MakeMinus(pc * Nc);
  ConvertMotionUpdates((float *) xms, &m_xv2s, &m_xba2s, &m_xbw2s);
  return true;
}


void LocalBundleAdjustor::EmbeddedPointIteration(const AlignedVector<Camera> &CsLF,
                                                 const AlignedVector<Rigid3D> &CsKF,
                                                 const std::vector<ubyte> &ucsKF,
                                                 const std::vector<ubyte> &uds,
                                                 std::vector<Depth::InverseGaussian> *ds) {
                                                   
  std::vector<int> &iKF2X = m_idxsTmp1, &iX2d = m_idxsTmp2;
  const int nKFs = static_cast<int>(m_KFs.size());
  iKF2X.assign(nKFs, -1);
  iX2d.resize(0);

  int Nd = 0;
  for (int iKF = 0; iKF < nKFs; ++iKF) {
    if (!(ucsKF[iKF] & LBA_FLAG_FRAME_UPDATE_DEPTH)) {
      continue;
    }
    const ubyte *_uds = uds.data() + m_iKF2d[iKF];
    const int iX = static_cast<int>(iX2d.size()), Nx = static_cast<int>(m_KFs[iKF].m_xs.size());
    iKF2X[iKF] = iX;
    iX2d.resize(iX + Nx, -1);
    int *ix2d = iX2d.data() + iX;
    for (int ix = 0; ix < Nx; ++ix) {
      if (_uds[ix] & LBA_FLAG_TRACK_UPDATE_DEPTH) {
        ix2d[ix] = Nd++;
      }
    }
  }

  int Nt = 0;
  m_idxsTmp3.resize(Nd + Nd + 1);
  int *Nzs = m_idxsTmp3.data(), *id2z = Nzs + Nd;
  memset(Nzs, 0, sizeof(int) * Nd);
  const int Nc = nKFs + static_cast<int>(m_LFs.size());
  for (int ic = 0; ic < Nc; ++ic) {
    const FRM::Frame *F = ic < nKFs ? (FRM::Frame *) &m_KFs[ic] : &m_LFs[ic - nKFs];
    const int NZ = static_cast<int>(F->m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = F->m_Zs[iZ];
      const int iX = iKF2X[Z.m_iKF];
      if (iX == -1) {
        continue;
      }
      bool t = false;
      const int *ix2d = iX2d.data() + iX;
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = F->m_zs[iz];
        const int id = ix2d[z.m_ix];
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
  for (int id = 0; id < Nd; ++id) {

    id2z[id + 1] = id2z[id] + Nzs[id];
  }
  m_zds.resize(id2z[Nd]);

  LA::Vector3f Rx;
  //LA::SymmetricMatrix2x2f W;
  Nt = 0;
  memset(Nzs, 0, sizeof(int) * Nd);
  for (int ic = 0; ic < Nc; ++ic) {
    const FRM::Frame *F = ic < nKFs ? (FRM::Frame *) &m_KFs[ic] : &m_LFs[ic - nKFs];
    const KeyFrame *KF = ic < nKFs ? (KeyFrame *) F : NULL;
    const LocalFrame *LF = KF ? NULL : (LocalFrame *) F;

    const Rigid3D &C = KF ? CsKF[ic] : CsLF[ic - nKFs].m_T;
    const int NZ = static_cast<int>(F->m_Zs.size());
    for (int iZ = 0; iZ < NZ; ++iZ) {
      const FRM::Measurement &Z = F->m_Zs[iZ];
      const int iX = iKF2X[Z.m_iKF];
      if (iX == -1) {
        continue;
      }
      const int *ix2d = iX2d.data() + iX;
      bool found = false;
      for (int iz = Z.m_iz1; iz < Z.m_iz2 && !found; ++iz) {
        found = ix2d[F->m_zs[iz].m_ix] != -1;
      }
      if (!found) {
        continue;
      }
      const Rigid3D T = C / CsKF[Z.m_iKF];
      LA::AlignedVector3f *t = m_t12s.Data() + Nt++;
      T.GetTranslation(*t);

      const KeyFrame &_KF = m_KFs[Z.m_iKF];
      for (int iz = Z.m_iz1; iz < Z.m_iz2; ++iz) {
        const FTR::Measurement &z = F->m_zs[iz];
        const int id = ix2d[z.m_ix];
        if (id == -1) {
          continue;
        }
        T.ApplyRotation(_KF.m_xs[z.m_ix].m_x, Rx);

        //z.m_W.GetScaled(KF ? KF->m_Azs[iz].m_wx : LF->m_Lzs[iz].m_wx, W);
        const LA::SymmetricMatrix2x2f &W = z.m_W;
        const int i = ++Nzs[id];
        m_zds[i].Set(*t, Rx, z.m_z, W);

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
    const KeyFrame &KF = m_KFs[iKF];
    const int Nx = static_cast<int>(KF.m_xs.size());
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
