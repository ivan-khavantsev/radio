package modem;

import audio.analysis.FFT;
import audio.analysis.FourierTransform;
import utils.Plot;

import java.awt.*;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;

public class Demodulator {

    private InputStream inputStream;
    FourierTransform fft = new FFT(Modulator.SAMPLES_PER_BIT, Modulator.SAMPLE_RATE);

    public Demodulator(InputStream inputStream) {
        this.inputStream = inputStream;
    }

    /*
        @int packetSize - size of data packet to demodulate
     */
    public synchronized byte[] demodulate(int packetSize) throws Exception {
        ByteBuffer dataBuffer = ByteBuffer.allocate(packetSize);
        this.sync();
        for (int i = 0; i < packetSize; i++) {
            byte b = this.readByte();
            dataBuffer.put(b);
        }
        return dataBuffer.array();
    }


    private void sync() throws IOException {
        /*
            Синхронизация последовательностью Бэкера
            sync = [+1, +1, +1, -1, -1, +1, -1];
         */
        byte[] syncBytes = new byte[Modulator.SYNC_BYTES.length];
        inputStream.read(syncBytes, 0, syncBytes.length);
        float[] audioFloats = this.getFloatSamples(syncBytes);


        while (true) {
            float[] rf = getRF(audioFloats);
            int max = max(rf); //????
            int symbol = (int)Math.floor(max/Modulator.SYNC_SAMPLES_PER_SYMBOL);
            if (rf[max]> 0.5f && symbol == 0) {
                break;
            } else {
                audioFloats = moveLeft(audioFloats, 3);
                byte[] tmp = new byte[6];
                inputStream.read(tmp, 0, tmp.length);
                float[] ftmp = getFloatSamples(tmp);
                System.arraycopy(ftmp, 0, audioFloats, audioFloats.length - 3, 3);
            }
        }
    }

    Plot plot = new Plot("Note A Spectrum", 1280, 512);
    public float[] getRF(float[] samples) {
        float[] rf = new float[samples.length]; // ??? - Значения взаимной корреляционной функции
        float[] sums = new float[ Modulator.SYNC_PACKET.length];
        float r = 0;
        for (int i = 0; i <  Modulator.SYNC_PACKET.length-1; i++) {
            float s = 0;
            for (int j = 0; j < Modulator.SYNC_PACKET.length-1; j++) {
                s = s + samples[i * Modulator.SYNC_SAMPLES_PER_SYMBOL + j];
            }
            sums[i] = s;
            r = r + (sums[i] * Modulator.SYNC_PACKET[i]) / Modulator.SYNC_SAMPLES_PER_SYMBOL;
        }

        rf[0] = r / Modulator.SYNC_PACKET.length;


        for (int i = 2; i < (samples.length - Modulator.SYNC_SAMPLES_PER_SYMBOL * Modulator.SYNC_PACKET.length); i++) {
            r = 0;
            for (int j = 1; j < Modulator.SYNC_PACKET.length; j++) {
                sums[j] = sums[j] - samples[i - 1 + (j - 1) * Modulator.SYNC_SAMPLES_PER_SYMBOL] + samples[i - 1 + j * Modulator.SYNC_SAMPLES_PER_SYMBOL];
                r = r + sums[j] * Modulator.SYNC_PACKET[j] / Modulator.SYNC_SAMPLES_PER_SYMBOL;
            }
            rf[i] = r / Modulator.SYNC_PACKET.length;
        }
        //plot.clear();
        plot.plot(samples, 1, Color.blue);
        plot.plot(rf, 1, Color.red);
        return rf;
    }

    public static float[] moveLeft(float[] array, int positions) {
        float[] temp = new float[array.length];
        System.arraycopy(array, positions, temp, 0, array.length - positions);
        return temp;
    }


    private byte readByte() throws IOException {
        int leftForByte = 8;
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < 8; i++) {
            byte[] buffer = new byte[Modulator.SAMPLES_PER_BIT * Modulator.BYTES_PER_SAMPLE];
            inputStream.read(buffer, 0, buffer.length);
            float[] audioFloats = this.getFloatSamples(buffer);
            int bit = this.getSamplesBit(audioFloats);

            sb.append(bit);

            leftForByte--;

            if (leftForByte < 1) {
                byte[] bval = new BigInteger(sb.toString(), 2).toByteArray();
                if (bval.length == 2) {
                    return bval[1];
                } else {
                    return bval[0];
                }
            }
        }
        return 0;
    }

    private synchronized int getSamplesBit(float[] audioFloats) {
        fft.forward(audioFloats);
        float[] spectrum = fft.getSpectrum();

        int bit;
        if (spectrum[0] > spectrum[1]) {
            bit = 0;
        } else {
            bit = 1;
        }
        return bit;
    }


    public static float[] getFloatSamples(byte[] bytes) {
        ShortBuffer sbuf = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN).asShortBuffer();
        short[] audioShorts = new short[sbuf.capacity()];
        sbuf.get(audioShorts);
        float[] audioFloats = new float[audioShorts.length];
        for (int i = 0; i < audioShorts.length; i++) {
            audioFloats[i] = ((float) audioShorts[i]) / 0x8000;
        }
        return audioFloats;
    }


    public static int max(float[] array) {
        int max = 0;
        for (int i = 0; i < array.length; i++) if (array[i] > array[max]) max = i;
        return max;
    }
}
