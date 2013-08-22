package modem;

import audio.analysis.FFT;
import audio.analysis.FourierTransform;
import utils.Utils;

import java.io.IOException;
import java.io.InputStream;

public class Demodulator {

    private InputStream inputStream;
    private FourierTransform fft = new FFT(Modulator.SAMPLES_PER_DATA_BIT, Modulator.SAMPLE_RATE);

    public Demodulator(InputStream inputStream) {
        this.inputStream = inputStream;
    }
    private boolean sync = false;
    public synchronized byte[] demodulate(int packetSize) throws Exception {
        byte[] withoutSyncBytes = new byte[Modulator.SYNC_READY_BYTES.length];
        if(!sync) {
            this.sync();
            sync = true;
        }else{
            inputStream.read(withoutSyncBytes);
        }

        byte[] data = this.readBytes(packetSize);
        return data;
    }

    public void resync(){
        sync = false;
    }


    private void sync() throws IOException {
        byte[] syncBytes = new byte[Modulator.SYNC_READY_BYTES.length];
        inputStream.read(syncBytes, 0, syncBytes.length);
        float[] audioFloats = Utils.bytesToFloats(syncBytes);

        while (true) {
            float maxAmplitude = audioFloats[Utils.max(audioFloats)];
            if (maxAmplitude > 0.3f && this.getSyncCorrelation(audioFloats) > maxAmplitude - 0.1f) {
                break;
            } else {
                audioFloats = Utils.moveLeft(audioFloats, 1);
                byte[] tempSampleBytes = new byte[Modulator.BYTES_PER_AUDIO_SAMPLE];
                inputStream.read(tempSampleBytes, 0, tempSampleBytes.length);
                float[] tempSample = Utils.bytesToFloats(tempSampleBytes);
                audioFloats[audioFloats.length - 1] = tempSample[0];
            }
        }
    }

    private byte[] readBytes(int count) throws IOException {
        byte[] result = new byte[count];
        byte[] buffer = new byte[Modulator.SAMPLES_PER_DATA_BIT * Modulator.BYTES_PER_AUDIO_SAMPLE];
        float[] audioFloats;
        for (int i = 0; i < count; i++) {
            byte b = 0;
            for (int j = 8; j > 0; j--) {
                inputStream.read(buffer, 0, buffer.length);
                audioFloats = Utils.bytesToFloats(buffer);
                int bit = this.samplesToBit(audioFloats);
                b = (byte) (b | bit << (j - 1));
            }
            result[i] = b;
        }
        return result;
    }

    private int samplesToBit1(float[] samples) {
        fft.forward(samples);

        float[] spectrum = fft.getSpectrum();

        int bit;
        if (spectrum[0] > spectrum[1]) {
            bit = 0;
        } else {
            bit = 1;
        }
        return bit;
    }

    private int samplesToBit(float[] samples) {

//        fft.forward(samples);
//        float[] spect = fft.getSpectrum();
//        fft.setBand(0,0);
//
//        for (int i = 2; i< spect.length;i++){
//            fft.setBand(i,0);
//        }
//
//        fft.inverse(samples);

        float corr = getCorrelation(samples,Modulator.ZERO_SAMPLES);
        if(corr > 0){
            return 0;
        }else{
            return 1;
        }
    }


    public static float getCorrelation(float[] signal, float[] etalon) {
        float signalSum = sum(signal, false);
        float etalonSum = sum(etalon, false);

        float signalXetalonSum = 0;
        for (int i = 0; i < signal.length; i++) {
            signalXetalonSum += signal[i] * etalon[i];
        }
        float top = (signal.length * signalXetalonSum) - (signalSum * etalonSum);
        float signalQuadSum = sum(signal, true);
        float etalonQuadSum = sum(etalon, true);

        float bottom1 = signal.length * signalQuadSum - signalSum*signalSum;
        float bottom2 = signal.length * etalonQuadSum - etalonSum*etalonSum;

        float bottom = (float)Math.sqrt(bottom1*bottom2);

        float result = top/bottom;
        return result;


    }

    public static float sum(float[] values, boolean quad) {
        float sum = 0;
        for (float v : values) {
            if (quad) {
                sum += v * v;
            } else {
                sum += v;
            }
        }
        return sum;
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
