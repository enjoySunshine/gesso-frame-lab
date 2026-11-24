/* Fused sRGB->linear conversion and integer box shrink.
 *
 * Linearising an 8-bit sRGB image to float scRGB and only then box-shrinking
 * it needs a full-resolution float image (4x the source memory). For
 * thumbnailing we can instead fold the sRGB->scRGB lookup into the box
 * average: every output pixel is the mean of the linearised values of its
 * block. The linearisation and the shrink therefore see the same pixels and
 * no full-resolution float buffer is ever allocated.
 *
 * The arithmetic mirrors running vips_sRGB2scRGB() followed by
 * vips_shrinkv()/vips_shrinkh(): the per-channel LUT is vips_v2Y_8[], the
 * block is accumulated in double and scaled by 1/n, then stored as float.
 */

/*

	This file is part of VIPS.

	VIPS is free software; you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
	02110-1301  USA

 */

/*

	These files are distributed with VIPS - http://www.vips.ecs.soton.ac.uk

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/
#include <glib/gi18n-lib.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include <vips/vips.h>
#include <vips/internal.h>

#include "../colour/pcolour.h"
#include "presample.h"

typedef struct _VipsLinearShrink {
	VipsResample parent_instance;

	int hshrink;
	int vshrink;

	/* Number of double accumulators in a per-sequence sum buffer.
	 */
	size_t n_sum;

} VipsLinearShrink;

typedef VipsResampleClass VipsLinearShrinkClass;

G_DEFINE_TYPE(VipsLinearShrink, vips_linear_shrink, VIPS_TYPE_RESAMPLE);

typedef struct {
	VipsRegion *ir;

	/* Block sums for the current chunk: double, bands * tile-width *
	 * chunk-height elements.
	 */
	double *sum;

} VipsLinearShrinkSequence;

static int
vips_linear_shrink_stop(void *vseq, void *a, void *b)
{
	VipsLinearShrinkSequence *seq = (VipsLinearShrinkSequence *) vseq;

	VIPS_UNREF(seq->ir);
	VIPS_FREE(seq->sum);
	VIPS_FREE(seq);

	return 0;
}

static void *
vips_linear_shrink_start(VipsImage *out, void *a, void *b)
{
	VipsImage *in = (VipsImage *) a;
	VipsLinearShrink *shrink = (VipsLinearShrink *) b;
	VipsLinearShrinkSequence *seq;

	if (!(seq = VIPS_NEW(NULL, VipsLinearShrinkSequence)))
		return NULL;

	seq->ir = vips_region_new(in);

	/* Number of double accumulators for the largest tile we might be
	 * asked to generate. In SMALLTILE mode a region is at most
	 * tile-width by fatstrip-height; size for the full output width to
	 * cover the wider FATSTRIP case too.
	 */
	shrink->n_sum =
		(size_t) out->Xsize * in->Bands * vips__fatstrip_height;
	seq->sum = VIPS_ARRAY(NULL, shrink->n_sum, double);
	if (!seq->sum) {
		vips_linear_shrink_stop(seq, a, b);
		return NULL;
	}

	return seq;
}

/* Generate one chunk of output. Each output pixel is the mean of the
 * linearised (v2Y LUT) values of its hshrink x vshrink input block.
 */
static int
vips_linear_shrink_gen(VipsRegion *out_region,
	void *vseq, void *a, void *b, gboolean *stop)
{
	VipsLinearShrink *shrink = (VipsLinearShrink *) b;
	VipsLinearShrinkSequence *seq =
		(VipsLinearShrinkSequence *) vseq;
	VipsImage *in = (VipsImage *) a;
	VipsRegion *ir = seq->ir;
	VipsRect *r = &out_region->valid;

	const int bands = in->Bands;
	const int hs = shrink->hshrink;
	const int vs = shrink->vshrink;
	const int dy = vips__fatstrip_height;

	vips_col_make_tables_RGB_8();

	for (int y = 0; y < r->height; y += dy) {
		int chunk_height = VIPS_MIN(dy, r->height - y);
		const int start = (r->top + y) * vs;
		const int end = (r->top + y + chunk_height) * vs;

		memset(seq->sum, 0,
			sizeof(double) * bands * r->width * chunk_height);

		VipsRect s;

		s.left = r->left * hs;
		s.top = start;
		s.width = r->width * hs;
		s.height = end - start;
		if (vips_region_prepare(ir, &s))
			return -1;

		/* Accumulate linearised values. Mirrors the float box shrink:
		 * vertical accumulation is done with a running double sum per
		 * output column, horizontal reduction inside each line.
		 */
		for (int yy = 0; yy < s.height; yy++) {
			int oy = (start + yy) / vs - (r->top + y);
			VipsPel *line = VIPS_REGION_ADDR(ir,
				s.left, s.top + yy);
			double *sp = seq->sum +
				(size_t) oy * r->width * bands;

			for (int x = 0; x < r->width; x++) {
				double *op = sp + (size_t) x * bands;

				for (int xx = 0; xx < hs; xx++) {
					VipsPel *pp =
						line + ((size_t) x * hs + xx) *
							bands;

					for (int bnd = 0; bnd < bands; bnd++)
						op[bnd] += vips_v2Y_8[pp[bnd]];
				}
			}
		}

		/* Average and write float output.
		 */
		const double inv = 1.0 / ((double) hs * vs);
		const int n = r->width * chunk_height * bands;
		float *q = (float *) VIPS_REGION_ADDR(out_region,
			r->left, r->top + y);

		for (int i = 0; i < n; i++)
			q[i] = (float) (seq->sum[i] * inv);
	}

	VIPS_COUNT_PIXELS(out_region, "vips_linear_shrink_gen");

	return 0;
}

static int
vips_linear_shrink_build(VipsObject *object)
{
	VipsObjectClass *class = VIPS_OBJECT_GET_CLASS(object);
	VipsResample *resample = VIPS_RESAMPLE(object);
	VipsLinearShrink *shrink = (VipsLinearShrink *) object;
	VipsImage **t = (VipsImage **) vips_object_local_array(object, 4);

	VipsImage *in;

	if (VIPS_OBJECT_CLASS(vips_linear_shrink_parent_class)->build(object))
		return -1;

	in = resample->in;

	if (shrink->hshrink < 1 ||
		shrink->vshrink < 1) {
		vips_error(class->nickname, "%s",
			_("shrink factors should be >= 1"));
		return -1;
	}

	if (in->BandFmt != VIPS_FORMAT_UCHAR) {
		vips_error(class->nickname, "%s",
			_("linear shrink needs uchar input"));
		return -1;
	}

	if (in->Bands != 3) {
		vips_error(class->nickname, "%s",
			_("linear shrink needs a 3-band sRGB image"));
		return -1;
	}

	/* No shrinking: linearise only, matching sRGB2scRGB.
	 */
	if (shrink->hshrink == 1 &&
		shrink->vshrink == 1)
		return vips_sRGB2scRGB(in, &t[0], NULL) ||
			vips_image_write(t[0], resample->out);

	/* Pad to a whole number of blocks so we never average half pixels,
	 * extending with edge copy just like shrinkv/shrinkh.
	 */
	if (in->Xsize % shrink->hshrink != 0 ||
		in->Ysize % shrink->vshrink != 0) {
		if (vips_embed(in, &t[1],
			0, 0,
			VIPS_ROUND_UP(in->Xsize, shrink->hshrink),
			VIPS_ROUND_UP(in->Ysize, shrink->vshrink),
			"extend", VIPS_EXTEND_COPY,
			NULL))
			return -1;
		in = t[1];
	}

	if (vips_image_pipelinev(resample->out,
		VIPS_DEMAND_STYLE_SMALLTILE, in, NULL))
		return -1;

	resample->out->BandFmt = VIPS_FORMAT_FLOAT;
	resample->out->Type = VIPS_INTERPRETATION_scRGB;

	/* Match the default (ceil == FALSE) sizing of shrinkh/shrinkv: round
	 * the divided dimension to nearest. The input is still padded to a
	 * whole number of blocks above so the edge block is averaged, never
	 * read out of bounds.
	 */
	resample->out->Xsize =
		VIPS_ROUND_UINT((double) resample->in->Xsize /
			shrink->hshrink);
	resample->out->Ysize =
		VIPS_ROUND_UINT((double) resample->in->Ysize /
			shrink->vshrink);

	if (resample->out->Xsize <= 0 ||
		resample->out->Ysize <= 0) {
		vips_error(class->nickname, "%s",
			_("image has shrunk to nothing"));
		return -1;
	}

	if (vips_image_generate(resample->out,
		vips_linear_shrink_start, vips_linear_shrink_gen,
		vips_linear_shrink_stop, in, shrink))
		return -1;

	return 0;
}

static void
vips_linear_shrink_class_init(VipsLinearShrinkClass *class)
{
	GObjectClass *gobject_class = G_OBJECT_CLASS(class);
	VipsObjectClass *vobject_class = VIPS_OBJECT_CLASS(class);
	VipsOperationClass *operation_class = VIPS_OPERATION_CLASS(class);

	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	vobject_class->nickname = "linearshrink";
	vobject_class->description =
		_("fused sRGB to linear conversion and box shrink");
	vobject_class->build = vips_linear_shrink_build;

	operation_class->flags |= VIPS_OPERATION_SEQUENTIAL;

	VIPS_ARG_INT(class, "hshrink", 3,
		_("Horizontal shrink"),
		_("Horizontal shrink factor"),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET(VipsLinearShrink, hshrink),
		1, VIPS_MAX_COORD, 1);

	VIPS_ARG_INT(class, "vshrink", 4,
		_("Vertical shrink"),
		_("Vertical shrink factor"),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET(VipsLinearShrink, vshrink),
		1, VIPS_MAX_COORD, 1);
}

static void
vips_linear_shrink_init(VipsLinearShrink *shrink)
{
	shrink->hshrink = 1;
	shrink->vshrink = 1;
}

/**
 * vips_linear_shrink:
 * @in: input 3-band uchar sRGB image
 * @out: (out): output float scRGB image
 * @hshrink: horizontal block shrink factor
 * @vshrink: vertical block shrink factor
 * @...: `NULL`-terminated list of optional named arguments
 *
 * Fuse an sRGB->linear (scRGB) conversion with an integer box shrink. Each
 * output sample is the mean, in linear light, of an @hshrink by @vshrink
 * block of input samples. The conversion and the average operate on the same
 * pixels, so the result is the linear-light shrink without ever holding a
 * full-resolution float image.
 *
 * Returns: 0 on success, -1 on error.
 */
int
vips_linear_shrink(VipsImage *in, VipsImage **out,
	int hshrink, int vshrink, ...)
{
	va_list ap;
	int result;

	va_start(ap, vshrink);
	result = vips_call_split("linearshrink", ap, in, out,
		hshrink, vshrink);
	va_end(ap);

	return result;
}
