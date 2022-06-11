import os

# import cv2
import imageio
from tqdm import tqdm

# source /vision2/u/yenyu/env-nerf/bin/activate

# These single object videos are already rendered because we didn't render frames in parallel.
# PRED: "/viscam/u/yenyu/osf-pytorch/logs/checkers_osf_{bunny,cup,bowl}_trans/renderonly_path_199999/video.mp4"
# GT: "/viscam/u/koven/osf/{bunny,cup,bowl}_relit_video"

OSF_DIRS = [
    "/viscam/u/yenyu/osf-pytorch/logs/checkers_comp_scene_1_0/renderonly_path_499999",
    "/viscam/u/yenyu/osf-pytorch/logs/checkers_comp_scene_1_nerf/renderonly_path_099999",
    "/viscam/u/yenyu/osf-pytorch/logs/checkers_osf_bunny_trans/renderonly_path_199999",
    "/viscam/u/yenyu/osf-pytorch/logs/checkers_osf_cup_trans/renderonly_path_199999",
    "/viscam/u/yenyu/osf-pytorch/logs/checkers_osf_bowl_trans/renderonly_path_199999",
]
PHYSG_DIRS = [
    "/viscam/u/mguo95/workspace/third_party/EricRyanChan/osf_siggraph_supp_vid/bunny",
    "/viscam/u/mguo95/workspace/third_party/EricRyanChan/osf_siggraph_supp_vid/cup",
    "/viscam/u/mguo95/workspace/third_party/EricRyanChan/osf_siggraph_supp_vid/bowl_trans",
]
GT_DIRS = [
    # "/viscam/u/koven/osf/comp_scene_3_new_const",
    "/viscam/u/koven/osf/bunny_relit_video",
    "/viscam/u/koven/osf/bunny_relit_video_3",  # larger
    "/viscam/u/koven/osf/bowl_relit_video_2",  # larger
    "/viscam/u/koven/osf/cup_relit_video",
    "/viscam/u/koven/osf/cup_relit_video_4",
    "/viscam/u/koven/osf/cup_relit_video_5",
    "/viscam/u/koven/osf/bowl_relit_video",
]


def main():

    # Generate OSF videos.
    # for img_dir in OSF_DIRS:
    #     mp4_path = os.path.join(img_dir, "gen_video.mp4")
    #     print(f"Generating video to {mp4_path}...")
    #     writer = imageio.get_writer(mp4_path, fps=8)
    #     for i in tqdm(range(40)):
    #         path = os.path.join(img_dir, f"{i:03}.png")
    #         writer.append_data(load_img(path))
    #     writer.close()
    
    # Generate PhySG videos.
    # for img_dir in PHYSG_DIRS:
    #     mp4_path = os.path.join(img_dir, "gen_video.mp4")
    #     print(f"Generating video to {mp4_path}...")
    #     writer = imageio.get_writer(mp4_path, fps=8)
    #     for i in tqdm(list(range(20, 40)) + list(range(20))):
    #         path = os.path.join(img_dir, f"sg_rgb_{i}.png")
    #         writer.append_data(load_img(path))
    #     writer.close()

    # Generate gt videos.
    for img_dir in GT_DIRS:
        mp4_path = os.path.join(img_dir, "gen_video.mp4")
        print(f"Generating video to {mp4_path}...")
        writer = imageio.get_writer(mp4_path, fps=8)
        paths = []
        for i in tqdm(range(0, 20)):  # 0 through 19
            path = os.path.join(img_dir, f"light_{i:03}.png")
            paths.append(path)
        for i in tqdm(range(1, 21)):  # 1 through 20
            path = os.path.join(img_dir, f"camera_{i:03}.png")
            paths.append(path)
        assert len(paths) == 40
        # fourcc = cv2.VideoWriter_fourcc(*'MPEG')
        # fps=8.0
        # fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
        # video = cv2.VideoWriter(mp4_path, fourcc, fps, (256, 256))
        # video = cv2.VideoWriter(mp4_path, fourcc, fps, (256, 256))
        # video = cv2.VideoWriter(mp4_path, 0, 1, (256, 256))
        for path in paths:
            if not os.path.exists(path):
                continue
            img = load_img(path)
            # video.write(img)

            writer.append_data(load_img(path))
        writer.close()
        # cv2.destroyAllWindows()
        # video.release()


def load_img(path):
    # BGR = cv2.imread(path, cv2.IMREAD_UNCHANGED)[:, :, :3]
    # RGB = cv2.cvtColor(BGR,cv2.COLOR_BGR2RGB)
    print(f"Loading path: {path}...")
    rgba = imageio.imread(path)[:, :, :3]
    # try:
    #     rgba = imageio.imread(path)
    # except:
    #     breakpoint()
    # rgb = rgba[:, :, 3]
    # imageio.imwrite(open(f"{path}_debug", "w"), rgb)
    # img = imageio.imread(path)[:, :, :3]
    # print(f"path: {path}")
    # print(f"rgb.shape: {rgb.shape}")
    return rgba


if __name__ == '__main__':
    main()
