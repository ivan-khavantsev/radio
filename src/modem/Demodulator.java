package modem;

import audio.analysis.FFT;
import audio.analysis.FourierTransform;

import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;

public class Demodulator {

    private InputStream inputStream;
    private FourierTransform fft = new FFT(Modulator.SAMPLES_PER_BIT, Modulator.SAMPLE_RATE);

    public Demodulator(InputStream inputStream) {
        this.inputStream = inputStream;
    }

    public synchronized byte[] demodulate(int packetSize) throws Exception {
        ByteBuffer dataBuffer = ByteBuffer.allocate(packetSize);
        this.sync1(packetSize);
        for (int i = 0; i < packetSize; i++) {
            byte b = this.readByte();
            dataBuffer.put(b);
        }
        return dataBuffer.array();
    }


    private void sync1(int packetSize) throws IOException {
        byte[] syncBytes = new byte[Modulator.SYNC_READY_BYTES.length];
        inputStream.read(syncBytes, 0, syncBytes.length);
        float[] audioFloats = this.getFloatSamples(syncBytes);


        while (true) {
            float maxAmplitude = audioFloats[max(audioFloats)];
            if (maxAmplitude > 0.3f && this.getSyncCorrelation(audioFloats) > maxAmplitude - 0.1f) {
                break;
            } else {
                audioFloats = moveLeft(audioFloats, 1);
                byte[] tmp = new byte[2];
                inputStream.read(tmp, 0, tmp.length);
                float[] newSamples = getFloatSamples(tmp);
                audioFloats[audioFloats.length - 1] = newSamples[0];
            }
        }
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
            byte[] buffer = new byte[Modulator.SAMPLES_PER_BIT * Modulator.BYTES_PER_AUDIO_SAMPLE];
            inputStream.read(buffer, 0, buffer.length);
            float[] audioFloats = this.getFloatSamples(buffer);
            int bit = this.getSamplesBit(audioFloats);

            sb.append(bit);

            leftForByte--;

            if (leftForByte < 1) {
                byte[] byteValue = new BigInteger(sb.toString(), 2).toByteArray();
                if (byteValue.length == 2) {
                    return byteValue[1];
                } else {
                    return byteValue[0];
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

    public static int max(float[] array, int offset, int length) {
        int max = offset;
        for (int i = offset; i < length; i++) if (array[i] > array[max]) max = i;
        return max;
    }

    public static int max(float[] array) {
        return max(array, 0, array.length);
    }


    public float getSyncCorrelation(float[] samples) {
        float[] sums = new float[Modulator.SYNC_SEQUENCE.length];
        float r = 0;
        for (int i = 0; i < Modulator.SYNC_SEQUENCE.length; i++) {
            float s = 0;
            for (int j = 0; j < Modulator.SYNC_SAMPLES_PER_SYMBOL; j++) {
                int cc = i * Modulator.SYNC_SAMPLES_PER_SYMBOL + j;
                s = s + samples[cc];
            }
            sums[i] = s;
            r = r + (sums[i] * Modulator.SYNC_SEQUENCE[i]) / Modulator.SYNC_SAMPLES_PER_SYMBOL;
        }

        float result = r / Modulator.SYNC_SEQUENCE.length;
        return result;
    }

}
