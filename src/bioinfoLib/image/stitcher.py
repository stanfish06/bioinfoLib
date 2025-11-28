# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

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
    img1_percentiles = np.percentile(img1[mask1], [1, 99])
    img2_percentiles = np.percentile(img2[mask2], [1, 99])
    img1_rescale = rescale_intensity(
        img1, in_range=(float(img1_percentiles[0]), float(img1_percentiles[1]))
    )
    img2_rescale = rescale_intensity(
        img2, in_range=(float(img2_percentiles[0]), float(img2_percentiles[1]))
    )
    detected_shift, _, _ = phase_cross_correlation(
        img1_rescale,
        img2_rescale,
        reference_mask=mask1,
        moving_mask=mask2,
        overlap_ratio=masked_overlap,
    )
    return detected_shift


# does not seem to work well for widefield
# confocal is a bit better
def linear_blend_2d(
    img1: np.ndarray,
    img2: np.ndarray,
    direction: Literal["up", "down", "left", "right"],
):
    sz1 = img1.shape
    sz2 = img2.shape
    assert sz1[0] == sz2[0] and sz1[1] == sz2[1], "input images must have same size"
    sz = sz1
    weight1 = np.ones_like(img1)
    match direction:
        case "up":
            weight1 = 1 - np.broadcast_to((np.arange(sz[0]) / (sz[0] - 1))[:, None], sz)
        case "down":
            weight1 = 1 - np.broadcast_to(
                (np.arange(sz[0])[::-1] / (sz[0] - 1))[:, None], sz
            )
        case "left":
            weight1 = 1 - np.broadcast_to((np.arange(sz[1]) / (sz[1] - 1))[None, :], sz)
        case "right":
            weight1 = 1 - np.broadcast_to(
                (np.arange(sz[1])[::-1] / (sz[1] - 1))[None, :], sz
            )
    return img1 * weight1 + img2 * (1 - weight1)


# runs quite slow for large image
def stitch_images_2d_one_direction(
    imgs,
    overlap: float,
    direction: Literal["up", "down", "left", "right"],
):
    if len(imgs) == 1:
        return imgs[0]

    # Calculate shifts between consecutive images
    shifts = []
    for i in trange(0, len(imgs) - 1, desc="Calculate shift"):
        shift = calculate_shift_2d_phase_correlation(
            imgs[i], imgs[i + 1], overlap=overlap, direction=direction
        )
        shifts.append(shift)

    shifts = np.array(shifts).astype(int)

    # Calculate canvas size and image positions
    positions = []
    canvas_size = [0, 0]

    match direction:
        case "up" | "down":
            # For vertical stitching
            offset_x = np.concatenate(([0], np.cumsum(shifts[:, 1])))
            offset_y = np.concatenate(([0], np.cumsum(shifts[:, 0])))

            # Calculate canvas width accounting for lateral shifts
            right_shift_max = np.max([abs(dx) for dx in offset_x if dx > 0] + [0])
            left_shift_max = np.max([abs(dx) for dx in offset_x if dx < 0] + [0])
            canvas_size[1] = imgs[0].shape[1] + right_shift_max + left_shift_max
            origin_x = left_shift_max

            # Calculate positions for each image
            for i in range(len(imgs)):
                pos_x = origin_x + offset_x[i]
                pos_y = offset_y[i] if direction == "down" else -offset_y[i]
                positions.append([pos_y, pos_x])

            # Adjust positions to be all positive and calculate canvas height
            min_y = min(pos[0] for pos in positions)
            for pos in positions:
                pos[0] -= min_y
            max_y = max(pos[0] + imgs[i].shape[0] for i, pos in enumerate(positions))
            canvas_size[0] = max_y

        case "left" | "right":
            # For horizontal stitching
            offset_x = np.concatenate(([0], np.cumsum(shifts[:, 1])))
            offset_y = np.concatenate(([0], np.cumsum(shifts[:, 0])))

            # Calculate canvas height accounting for vertical shifts
            up_shift_max = np.max([abs(dy) for dy in offset_y if dy < 0] + [0])
            down_shift_max = np.max([abs(dy) for dy in offset_y if dy > 0] + [0])
            canvas_size[0] = imgs[0].shape[0] + up_shift_max + down_shift_max
            origin_y = up_shift_max

            # Calculate positions for each image
            for i in range(len(imgs)):
                pos_x = offset_x[i] if direction == "right" else -offset_x[i]
                pos_y = origin_y + offset_y[i]
                positions.append([pos_y, pos_x])

            # Adjust positions to be all positive and calculate canvas width
            min_x = min(pos[1] for pos in positions)
            for pos in positions:
                pos[1] -= min_x
            max_x = max(pos[1] + imgs[i].shape[1] for i, pos in enumerate(positions))
            canvas_size[1] = max_x

    # Create canvas and place images
    stitched_img = np.zeros(canvas_size, dtype=imgs[0].dtype)

    # Place first image
    pos = positions[0]
    stitched_img[
        pos[0] : pos[0] + imgs[0].shape[0], pos[1] : pos[1] + imgs[0].shape[1]
    ] = imgs[0]

    # Place subsequent images with blending
    for i in range(1, len(imgs)):
        pos = positions[i]
        img_h, img_w = imgs[i].shape[:2]

        # Get the region where this image will be placed
        y1, y2 = pos[0], pos[0] + img_h
        x1, x2 = pos[1], pos[1] + img_w

        # Get current canvas region
        canvas_region = stitched_img[y1:y2, x1:x2]

        # Find overlap region (where canvas has non-zero values)
        overlap_mask = canvas_region > 0

        if np.any(overlap_mask):
            # Create images of same size for blending
            img1_blend = np.zeros_like(imgs[i])
            img2_blend = np.zeros_like(imgs[i])

            img1_blend[overlap_mask] = canvas_region[overlap_mask]
            img2_blend = imgs[i].copy()

            # Create a mask for the overlap region
            overlap_region = np.zeros_like(imgs[i])
            overlap_region[overlap_mask] = 1

            # Blend only in overlap region
            blended = linear_blend_2d(img1_blend, img2_blend, direction)

            # Place the image: non-overlapping parts directly, overlapping parts blended
            result_img = imgs[i].copy()
            result_img[overlap_mask] = blended[overlap_mask]

            stitched_img[y1:y2, x1:x2] = result_img
        else:
            # No overlap, place directly
            stitched_img[y1:y2, x1:x2] = imgs[i]

    return stitched_img
