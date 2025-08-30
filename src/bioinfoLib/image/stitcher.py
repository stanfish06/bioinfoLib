from typing import Literal

import numpy as np
from skimage.exposure import rescale_intensity
from skimage.registration import phase_cross_correlation
from tqdm import trange


def calculate_shift_2d_phase_correlation(
    img1,
    img2,
    overlap: float,
    direction: Literal["up", "down", "left", "right"],
    masked_overlap: float = 0.5,
):
    # create masks based on direction
    # assume img1 and img2 having same size
    mask1 = np.zeros_like(img1)
    mask2 = np.zeros_like(img1)
    size_y, size_x = img1.shape
    match direction:
        case "up":
            mask1[0 : int(np.floor(size_y * overlap)), :] = True
            mask2[int(np.ceil(size_y * (1 - overlap))) : size_y, :] = True
        case "down":
            mask1[int(np.ceil(size_y * (1 - overlap))) : size_y, :] = True
            mask2[0 : int(np.floor(size_y * overlap)), :] = True
        case "left":
            mask1[:, 0 : int(np.floor(size_x * overlap))] = True
            mask2[:, int(np.ceil(size_x * (1 - overlap))) : size_x] = True
        case "right":
            mask1[:, int(np.ceil(size_x * (1 - overlap))) : size_x] = True
            mask2[:, 0 : int(np.floor(size_x * overlap))] = True
    # adjust intensity based on overlapping region
    img1_rescale = rescale_intensity(img1, tuple(np.percentile(img1[mask1], [1, 99])))
    img2_rescale = rescale_intensity(img2, tuple(np.percentile(img2[mask2], [1, 99])))
    detected_shift, _, _ = phase_cross_correlation(
        img1_rescale,
        img2_rescale,
        reference_mask=mask1,
        moving_mask=mask2,
        overlap_ratio=masked_overlap,
    )
    return detected_shift


def stitch_images_2d_one_direction(
    imgs,
    overlap: float,
    direction: Literal["up", "down", "left", "right"],
):
    shifts = []
    canvas_size = [0, 0]
    match direction:
        case "up" | "down":
            canvas_size[1] = imgs[0].shape[1]
            # add last img's height
            canvas_size[0] = imgs[-1].shape[0]
        case "left" | "right":
            canvas_size[0] = imgs[0].shape[0]
            # add last img's width
            canvas_size[1] = imgs[-1].shape[1]
    for i in trange(0, len(imgs) - 1, desc="Calculate shift"):
        shift = calculate_shift_2d_phase_correlation(
            imgs[i], imgs[i + 1], overlap=overlap, direction=direction
        )
        shifts.append(shift)
        match direction:
            case "up" | "down":
                canvas_size[0] += abs(shift[0])
            case "left" | "right":
                canvas_size[1] += abs(shift[1])

    origin_xy = [0, 0]
    shifts = np.array(shifts).astype(int)
    match direction:
        case "up" | "down":
            offset_x = np.cumsum(shifts[:, 1]).astype(int)
            right_shift_max = np.max([abs(dx) for dx in offset_x if dx > 0] + [0])
            left_shift_max = np.max([abs(dx) for dx in offset_x if dx < 0] + [0])
            canvas_size[1] += right_shift_max + left_shift_max
            origin_xy[0] = left_shift_max
            tile_pos = np.cumsum(shifts[:, 0]).astype(int)
        case "left" | "right":
            offset_y = np.cumsum(shifts[:, 0]).astype(int)
            up_shift_max = np.max([abs(dy) for dy in offset_y if dy < 0] + [0])
            down_shift_max = np.max([abs(dy) for dy in offset_y if dy > 0] + [0])
            canvas_size[0] += up_shift_max + down_shift_max
            origin_xy[1] = up_shift_max
            tile_pos = np.cumsum(shifts[:, 1]).astype(int)
    stitched_img = np.zeros(np.array(canvas_size).astype(int), dtype=imgs[0].dtype)
    # TODO: blending in the overlay region. Ignore for now
    for i in range(len(imgs)):
        match direction:
            case "up":
                if i == 0:
                    stitched_img[
                        tile_pos[i] :,
                        (origin_xy[0]) : (origin_xy[0] + imgs.shape[1]),
                    ] = imgs[i][shifts[i][0] :, :]
                elif i < len(imgs):
                    stitched_img[
                        tile_pos[i] : tile_pos[i - 1],
                        (origin_xy[0] + offset_x[i - 1]) : (
                            origin_xy[0] + offset_x[i - 1][1] + imgs.shape[1]
                        ),
                    ] = imgs[i][shifts[i][0] :, :]
                else:
                    stitched_img[
                        0 : tile_pos[i - 1],
                        (origin_xy[0] + offset_x[i - 1]) : (
                            origin_xy[0] + offset_x[i - 1] + imgs.shape[1]
                        ),
                    ] = imgs[i]
            case "down":
                if i == 0:
                    stitched_img[
                        0 : tile_pos[i],
                        (origin_xy[0]) : (origin_xy[0] + imgs.shape[1]),
                    ] = imgs[i][0 : shifts[i][0], :]
                elif i < len(imgs) - 1:
                    stitched_img[
                        tile_pos[i - 1] : tile_pos[i],
                        (origin_xy[0] + offset_x[i - 1]) : (
                            origin_xy[0] + offset_x[i - 1] + imgs.shape[1]
                        ),
                    ] = imgs[i][0 : shifts[i][0], :]
                else:
                    stitched_img[
                        tile_pos[i - 1] :,
                        (origin_xy[0] + offset_x[i - 1]) : (
                            origin_xy[0] + offset_x[i - 1] + imgs.shape[1]
                        ),
                    ] = imgs[i]
            case "left":
                if i == 0:
                    stitched_img[
                        (origin_xy[1]) : (origin_xy[1] + imgs[i].shape[0]),
                        tile_pos[i] :,
                    ] = imgs[i][:, shifts[i][1] :]
                elif i < len(imgs) - 1:
                    stitched_img[
                        (origin_xy[1] + offset_y[i - 1]) : (
                            origin_xy[1] + offset_y[i - 1] + imgs[i].shape[0]
                        ),
                        tile_pos[i] : tile_pos[i - 1],
                    ] = imgs[i][:, shifts[i][1] :]
                else:
                    stitched_img[
                        (origin_xy[1] + offset_y[i - 1]) : (
                            origin_xy[1] + offset_y[i - 1] + imgs[i].shape[0]
                        ),
                        0 : tile_pos[i - 1],
                    ] = imgs[i]
            case "right":
                if i == 0:
                    stitched_img[
                        (origin_xy[1]) : (origin_xy[1] + imgs[i].shape[0]),
                        0 : tile_pos[i],
                    ] = imgs[i][:, 0 : shifts[i][1]]
                elif i < len(imgs) - 1:
                    stitched_img[
                        (origin_xy[1] + offset_y[i - 1]) : (
                            origin_xy[1] + offset_y[i - 1] + imgs[i].shape[0]
                        ),
                        tile_pos[i - 1] : tile_pos[i],
                    ] = imgs[i][:, 0 : shifts[i][1]]
                else:
                    stitched_img[
                        (origin_xy[1] + offset_y[i - 1]) : (
                            origin_xy[1] + offset_y[i - 1] + imgs[i].shape[0]
                        ),
                        tile_pos[i - 1] :,
                    ] = imgs[i]
    return stitched_img
