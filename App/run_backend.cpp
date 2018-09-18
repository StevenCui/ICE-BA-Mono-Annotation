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


// solver.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "InputSequence.h"
#include "Parameter.h"
#include "Timer.h"
#include "MultiThread.h"

#ifdef IBA_WITH_CVD
#include "ViewerIBA.h"
#endif

#include <gtest/gtest.h>
#ifndef WIN32
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <chrono>
#include <thread>
#elif defined IBA_WITH_CVD
#include "Z:/CVD/lib/link.h"
#include "Z:/Graphics/glut/lib/link.h"
#include "Z:/Graphics/glew/lib/link.h"
#include "Z:/GoogleTest/lib/link.h"
#endif
#include <string>
#include <vector>

bool RunSolver(const Configurator &cfgor, const InputSequence &IS, const std::string param = "",
               IBA::Time *tLBA = NULL, IBA::Time *tGBA = NULL
             ) {
  const int nFrms = IS.Size();
  if (nFrms == 0) {
    return false;
  }

  std::vector<int> seeds;
  const int seed = cfgor.GetArgument("random_seed", 0);
  if (seed != 0) {
    srand(seed);
    seeds.resize(nFrms);
    for (int iFrm = 0; iFrm < nFrms; ++iFrm) {
      seeds[iFrm] = UT::Random<int>();
    }
  }

  IBA::Solver solver;
  const int serial = (cfgor.GetArgument("serial_lba", IBA_SERIAL_NONE) & 255) |
                    ((cfgor.GetArgument("serial_gba", IBA_SERIAL_NONE) & 255) << 8);
  const int verbose = (cfgor.GetArgument("verbose_lba", IBA_VERBOSE_NONE) & 255) |
                     ((cfgor.GetArgument("verbose_gba", IBA_VERBOSE_NONE) & 255) << 8);
  const int debug = (cfgor.GetArgument("debug_lba", IBA_DEBUG_NONE) & 255) |
                   ((cfgor.GetArgument("debug_gba", IBA_DEBUG_NONE) & 255) << 8);
  const int history = (cfgor.GetArgument("history_lba", IBA_HISTORY_NONE) & 255) |
                     ((cfgor.GetArgument("history_gba", IBA_HISTORY_NONE) & 255) << 8);
  solver.Create(IS.m_K, serial, verbose, debug, history, param, IS.m_dir

              );
#ifndef WIN32
    if (cfgor.GetArgument("serial_lba", IBA_SERIAL_NONE) > 0) {
      LOG(INFO) << "LBA is working serial";
    } else {
      LOG(INFO) << "LBA is working parallel";
    }
    if (cfgor.GetArgument("serial_gba", IBA_SERIAL_NONE) > 0) {
      LOG(INFO) << "GBA is working serial";
    } else {
      LOG(INFO) << "GBA is working parallel";
    }
#endif
  solver.Start();
#ifdef IBA_WITH_CVD
  ViewerIBA viewer;
  const int visualize = cfgor.GetArgument("visualize", 1);
  viewer.Create(&solver, cfgor.GetArgument("screen_file"),
                         cfgor.GetArgument("screen_combine", 0),
                         cfgor.GetArgument("save_file"),
                         cfgor.GetArgument("save_frame", INT_MAX), visualize != 0);
  const int pause = cfgor.GetArgument("pause", 0);
  const int iFrmStart = viewer.Start(cfgor.GetArgument("view_file_load"), pause >= 2);
#else
  const int iFrmStart = -1;
#endif


  IBA::CurrentFrame CF;
  IBA::KeyFrame KF;
  IBA::RelativeConstraint Z;
  std::vector<int> iFrms;
  std::vector<IBA::CameraPose> Cs;
  IBA::SlidingWindow SW;
  std::vector<int> idxs;

  // IBA::Error eLBA, eGBA;
  const int print = cfgor.GetArgument("print_camera", 0);

  const float dkfRatio = cfgor.GetArgument("delete_keyframe_ratio", 0.0f);
  const float dmpRatioFrq = cfgor.GetArgument("delete_map_point_ratio_frequency", 0.0f);
  const float dmpRatioNum = cfgor.GetArgument("delete_map_point_ratio_number", 0.0f);
  //UT::PrintStart(UT::FileNameAppendSuffix(UT::FileNameReplaceDirectory(
  //               cfgor.GetArgument("print_file"), ".", IS.m_dir)));
  UT::PrintStart(UT::FileNameReplaceDirectory(cfgor.GetArgument("print_file"),
                 ".", IS.m_dir));
  for (int iFrm = iFrmStart + 1, iFrmLast = -1; iFrm < nFrms; ++iFrm) {
    const double t1 = Timer::GetTime();
    if (!seeds.empty()) {
      srand(seeds[iFrm]);
    }

    IS.LoadCurrentFrame(iFrm, &CF, &KF, &solver);

    solver.PushCurrentFrame(CF, KF.iFrm == -1 ? NULL : &KF);

    if (IS.LoadRelativeConstraint(iFrm, &Z)) {

      S.MakeDiagonal(s2p, s2r);
      S.Get(Z.S.S);
      solver.PushRelativeConstraint(Z);
    }
    if (IS.LoadDeleteKeyFrames(iFrm, &iFrms)) {
      const int N = static_cast<int>(iFrms.size());
      for (int i = 0; i < N; ++i) {
        const int jFrm = iFrms[i];
        solver.DeleteKeyFrame(jFrm);
#ifdef IBA_WITH_CVD
        viewer.DeleteKeyFrame(jFrm);
#endif
      }
    }
    if (IS.LoadDeleteMapPoints(iFrm, &idxs)) {
      solver.DeleteMapPoints(idxs);
    }
    //UT::Print("[%d] %f\n", iFrm, Timer::GetTime() - t1);
    if (IS.LoadUpdateCameras(iFrm, &iFrms, &Cs)) {
      solver.UpdateCameras(iFrms, Cs);
    }
    const bool scc = solver.GetSlidingWindow(&SW);
    if (scc) {
      //IBA::PrintPoints(SW.Xs);
      if (print > 0) {
        std::vector<int>::const_iterator i;
        for (i = std::upper_bound(SW.iFrms.begin(), SW.iFrms.end(), iFrmLast);
             i != SW.iFrms.end(); ++i) {
          iFrmLast = *i;
          const int j = static_cast<int>(i - SW.iFrms.begin());
          IBA::PrintCameraPose(iFrmLast, SW.CsLF[j].C, print > 1);
        }
      }
    }

    if (dkfRatio > 0.0f && UT::Random<float>() < dkfRatio) {
      //const int iFrm0 = SW.iFrms.front();
      const int iFrm0 = iFrm - LBA_MAX_LOCAL_FRAMES + 1;
      const int nKFs = solver.SearchKeyFrame(iFrm0);
      if (nKFs > 1) {
        const int i = UT::Random<int>(nKFs), jFrm = solver.GetKeyFrameIndex(i);
        solver.DeleteKeyFrame(jFrm);
#ifdef IBA_WITH_CVD
        viewer.DeleteKeyFrame(jFrm);
#endif
      }
    }
    if (dmpRatioFrq > 0.0f && UT::Random<float>() < dmpRatioFrq) {
      solver.GetMapPointIndexes(&idxs);
      int i, j;
      const int N = static_cast<int>(idxs.size());
      for (i = j = 0; i < N; ++i) {
        if (UT::Random<float>() < dmpRatioNum) {
          idxs[j++] = idxs[i];
        }
      }
      idxs.resize(j);
      solver.DeleteMapPoints(idxs);
    }
    if (print == 0) {
      printf("\r%d / %d = %f%%", iFrm + 1, nFrms, UT::Percentage(iFrm + 1, nFrms));
    }
#ifdef IBA_WITH_CVD
    if (scc) {
      if (!viewer.Run(visualize >= 2, true, iFrm)) {
        break;
      }
    }
#endif
    const double t2 = Timer::GetTime(), dt = t2 - t1;

    if (dt < IS.m_dt) {
      const int ms = static_cast<int>((IS.m_dt - dt) * 1000.0 + 0.5);
#ifdef WIN32
      Sleep(ms);  // Windows Sleep (ms)
#else
      std::this_thread::sleep_for(std::chrono::milliseconds(ms));
#endif
    }
  }
  if (print == 0 || print == 1) {
    printf("\n");
  }
  solver.Stop();

  UT::Print("%.2f\t%.2f * %d\n", tLBA->t / tLBA->n, tGBA->t / tGBA->n, tGBA->n);


  const bool append = cfgor.GetArgument("output_append", 1) != 0;
  const bool poseOnly = cfgor.GetArgument("output_camera_pose_only", 1) != 0;
  solver.SaveCamerasLBA(cfgor.GetArgument("output_camera_file_lba"), append, poseOnly);
  solver.SaveCamerasGBA(cfgor.GetArgument("output_camera_file_gba"), append, poseOnly);
  solver.SaveCostsLBA(cfgor.GetArgument("output_cost_file_lba_a"), append, 0);
  solver.SaveCostsGBA(cfgor.GetArgument("output_cost_file_gba_a"), append, 0);
  solver.SaveCostsLBA(cfgor.GetArgument("output_cost_file_lba_b"), append, 1);
  solver.SaveCostsGBA(cfgor.GetArgument("output_cost_file_gba_b"), append, 1);
  solver.SaveCostsLBA(cfgor.GetArgument("output_cost_file_lba_p"), append, 2);
  solver.SaveCostsGBA(cfgor.GetArgument("output_cost_file_gba_p"), append, 2);

  solver.SaveResidualsLBA(cfgor.GetArgument("output_residual_file_lba"), append);
  solver.SaveResidualsGBA(cfgor.GetArgument("output_residual_file_gba"), append);

  solver.SaveMarginalizations(cfgor.GetArgument("output_marg_file_lp"), append, 0);
  solver.SaveMarginalizations(cfgor.GetArgument("output_marg_file_em"), append, 1);

  solver.SavePointsLBA(cfgor.GetArgument("output_point_file_lba"), append);
  solver.SavePointsGBA(cfgor.GetArgument("output_point_file_gba"), append);
  if (tLBA) {
    solver.GetTimeLBA(tLBA);
  }
  if (tGBA) {
    solver.GetTimeGBA(tGBA);
  }
  solver.SaveTimesLBA(cfgor.GetArgument("output_time_file_lba"), append);
  solver.SaveTimesGBA(cfgor.GetArgument("output_time_file_gba"), append);

  if (tLBA && tGBA) {

    UT::Print("%.2f\t%.2f * %d\n", tLBA->t / tLBA->n, tGBA->t / tGBA->n, tGBA->n);

  }
  UT::PrintStop();
#ifdef IBA_WITH_CVD
  viewer.Stop(cfgor.GetArgument("view_file_save"), visualize == 1 || pause >= 1);
#endif

  solver.Destroy();
  return true;
}

bool RunConfig(const std::string &fileName) {
  Configurator cfgor(fileName);
  IBA::Time tLBA, tGBA;
  const std::string dir = cfgor.GetArgument("input_directory_all");
  if (dir == "") {
    const InputSequence IS(cfgor);
    return RunSolver(cfgor, IS, fileName, &tLBA, &tGBA);
  } else {
    IBA::Time StLBA, StGBA;
    StLBA.t = StGBA.t = 0.0f;
    StLBA.n = StGBA.n = 0;
    // StLBA.nd = StGBA.nd = 0;
    std::vector<std::string> seqName;

    const std::string listName = dir + cfgor.GetArgument("input_sequence_list", "list.txt");
    FILE *fp = fopen(listName.c_str(), "r");
    if (!fp) {
      return false;
    }
    char line[UT_STRING_WIDTH_MAX];
    while (fgets(line, UT_STRING_WIDTH_MAX, fp)) {
      if (line[0] == '#') {
        continue;
      }
      const int len = static_cast<int>(strlen(line));
      if (line[len - 1] == 10) {
        line[len - 1] = 0;
      }
      seqName.push_back(line);
    }
    fclose(fp);
    UT::PrintLoaded(listName);

    cfgor.SetArgument("save_file", "");
    cfgor.SetArgument("visualize", 0);
    cfgor.SetArgument("pause", 0);
    cfgor.SetArgument("verbose_lba", 0);
    cfgor.SetArgument("verbose_gba", 0);
    const int N = static_cast<int>(seqName.size());
    for (int i = 0; i < N; ++i) {


      UT::Print("%s\t", seqName[i].c_str());

      cfgor.SetArgument("input_directory", dir + seqName[i] + "/");
// #ifdef CFG_GROUND_TRUTH
//    cfgor.SetArgument("input_ground_truth_time_first", tFirstGT[i]);
// #endif
      const InputSequence IS(cfgor);

      if (!RunSolver(cfgor, IS, fileName, &tLBA, &tGBA)) {

        return false;
      }

      StLBA.t += tLBA.t;  StLBA.n += tLBA.n;
      StGBA.t += tGBA.t;  StGBA.n += tGBA.n;
      // StLBA.nd += tLBA.nd;
      // StGBA.nd += tGBA.nd;
    }

    UT::Print("Total\t\t%.2f\t%.2f * %d\n", StLBA.t / StLBA.n, StGBA.t / StGBA.n, StGBA.n);

    return true;
  }
}

#ifdef WIN32
int _tmain(int argc, _TCHAR* argv[]) {
#else
int main(int argc, char** argv) {
#endif  // _WIND32
// #ifdef CFG_DEBUG
#ifndef WIN32
  google::InitGoogleLogging(argv[0]);
  // FLAGS_stderrthreshold = 0;
#endif
#ifndef WIN32
  google::ParseCommandLineFlags(&argc, &argv, false);
#endif
  // TODO(mingyu): Upgrade to use gflag once it's ready on Windows
  if (argc == 2) {
    // Playback according to the input config file
    const bool scc = RunConfig(std::string(argv[1]));
    if (scc) {
      return 0;
    } else {
      return 1;
    }
  } else {
    std::cout << " Incorrect input arguments\n";
    return 0;
  }
}
