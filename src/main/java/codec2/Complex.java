/*
 * Copyright (C) 1990-2013 David Rowe, VK5DGR
 *
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 2.1, as published
 * by the Free Software Foundation. This program is distributed in the hope that
 * it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */
package codec2;

/**
 * Complex using 32-bit floats. This class consumes a lot of time and memory.
 *
 * <p>
 * When declaring Complex arrays, you must also allocate the storage. For
 * example, after declaring:
 *
 * <p>
 * private Complex[] input = new Complex[FRAME_LENGTH];
 *
 * <p>
 * You must then create the objects:
 *
 * <p>
 * for (int i = 0; i < FRAME_LENGTH; i++) {<br> input[i] = new Complex();<br> }
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Complex {

    private final float re;
    private final float im;

    public Complex() {
        this.re = 0.0F;
        this.im = 0.0F;
    }

    public Complex(float real, float imag) {
        this.re = real;
        this.im = imag;
    }

    public float getReal() {
        return this.re;
    }

    public float getImaginary() {
        return this.im;
    }

    public Complex add(Complex b) {
        return new Complex(this.re + b.re, this.im + b.im);
    }

    public Complex minus(Complex b) {
        return new Complex(this.re - b.re, this.im - b.im);
    }

    public Complex times(Complex b) {
        return new Complex(this.re * b.re - this.im * b.im, this.re * b.im + this.im * b.re);
    }

    public Complex times(float alpha) {
        return new Complex(this.re * alpha, this.im * alpha);
    }

    public Complex conjugate() {
        return new Complex(this.re, -this.im);
    }

    public Complex divide(float val) {
        return new Complex(this.re / val, this.im / val);
    }

    public Complex cneg() {
        return new Complex(-this.re, -this.im);
    }

    public float csqr() {
        return ((this.re * this.re) + (this.im * this.im));
    }

    public float cabsolute() {
        return (float) Math.sqrt(csqr());
    }
}
