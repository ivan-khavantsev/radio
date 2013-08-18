import audio.AudioChannel;
import crypt.AES;
import modem.Demodulator;
import org.xiph.speex.SpeexDecoder;

import javax.sound.sampled.*;

public class Receiver implements Runnable {
    static SpeexDecoder speexDecoder = new SpeexDecoder();
    static SourceDataLine line2 = null;

    @Override
    public void run() {
        try {

            AudioFormat audioFormat2 = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,
                    8000.0F, 16, 1, 2, 8000.0F, false);
            Mixer.Info[] mixerInfo = AudioSystem.getMixerInfo();
            Mixer mixer = AudioSystem.getMixer(mixerInfo[2]);
            line2 = (SourceDataLine) mixer.getLine(new DataLine.Info(SourceDataLine.class, audioFormat2));
            line2.open();
            line2.start();


            speexDecoder.init(0, 8000, 1, true);
            int packetSize = 16;
            Demodulator demodulator = new Demodulator(new AudioInputStream(AudioChannel.getInLine()));


            AES aes = new AES("abc");
            while (true) {

                try{
                    byte[] encodedArray = demodulator.demodulate(packetSize);
                   // System.out.println(new String(encodedArray));

                    byte[] decrypted = aes.decrypt(encodedArray);

                    speexDecoder.processData(decrypted, 0, decrypted.length/2);
                    byte[] decodedData = new byte[speexDecoder.getProcessedDataByteSize()];
                    speexDecoder.getProcessedData(decodedData, 0);
                    line2.write(decodedData, 0, decodedData.length);

                    speexDecoder.processData(decrypted, decrypted.length/2, decrypted.length/2);
                    decodedData = new byte[speexDecoder.getProcessedDataByteSize()];
                    speexDecoder.getProcessedData(decodedData, 0);
                    line2.write(decodedData, 0, decodedData.length);




                }catch (Throwable t){
                    System.out.println(t);
                }

            }


        } catch (Throwable t) {
            System.out.println(t);
        }
    }
}
