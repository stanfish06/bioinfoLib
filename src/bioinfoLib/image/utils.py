# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import bioformats
import numpy as np
import tifffile
from tqdm import trange


def read_img_stack(path):
    meta = bioformats.get_omexml_metadata(path)
    meta = bioformats.omexml.OMEXML(meta).image(0)
    size_t = meta.Pixels.get_SizeT()
    size_z = meta.Pixels.get_SizeZ()
    size_y = meta.Pixels.get_SizeY()
    size_x = meta.Pixels.get_SizeX()
    size_c = meta.Pixels.get_SizeC()

    # not sure how to get time interval, leave it for now
    pixel_size_z = meta.Pixels.get_PhysicalSizeZ()
    pixel_size_y = meta.Pixels.get_PhysicalSizeY()
    pixel_size_x = meta.Pixels.get_PhysicalSizeX()

    unit_xy = meta.Pixels.get_PhysicalSizeXUnit()
    unit_z = meta.Pixels.get_PhysicalSizeZUnit()

    stack = np.zeros((size_t, size_z, size_c, size_y, size_x), dtype=np.uint16)
    with bioformats.ImageReader(path) as reader:
        for t in trange(size_t, desc="Read T-frames", position=0):
            for z in trange(size_z, desc=f"Read Z-slices (t = {t + 1})", leave=False):
                for c in range(size_c):
                    img = reader.read(t=t, z=z, c=c)
                    if img.dtype == np.float32 or img.dtype == np.float64:
                        img = (img * 65535).astype(np.uint16)
                    stack[t, z, c, ...] = img

    meta = {
        "pixel_size_x": pixel_size_x,
        "pixel_size_y": pixel_size_y,
        "pixel_size_z": pixel_size_z,
        "unit_xy": unit_xy,
        "unit_z": unit_z,
    }
    return stack, meta


def save_img_stack(path: str, img, meta):
    resolution_x = 1.0 / meta["pixel_size_x"] if meta["pixel_size_x"] else None
    resolution_y = 1.0 / meta["pixel_size_y"] if meta["pixel_size_y"] else None
    tifffile.imwrite(
        path,
        img,
        imagej=True,
        resolution=(resolution_x, resolution_y),
        resolutionunit=f"1/{meta['unit']}",
        metadata={
            "axes": "TZCYX",
            "spacing": meta["pixel_size_z"],
            "unit": meta["unit"],
        },
    )
