import audio.AudioChannel;
import crypt.AES;
import modem.Modulator;
import org.xiph.speex.SpeexDecoder;
import org.xiph.speex.SpeexEncoder;

import javax.sound.sampled.*;

public class Sender implements Runnable {
    AudioFormat channelAudioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, 44100, 16, 1, 2, 44100, false);
    AudioFormat micAudioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, 8000.0F, 16, 1, 2, 8000.0F, false);

    @Override
    public void run() {
        try {


            Mixer.Info[] mixerInfo = AudioSystem.getMixerInfo();
            Mixer mixer = AudioSystem.getMixer(mixerInfo[7]);
            final TargetDataLine line = (TargetDataLine) mixer.getLine(new DataLine.Info(TargetDataLine.class, micAudioFormat));
            line.open();
            line.start();


            SpeexEncoder speexEncoder = new SpeexEncoder();
            speexEncoder.init(0, 0, 8000, 1);
            int packetSize = speexEncoder.getFrameSize() * 2;

            byte[] inBuffer = new byte[packetSize];

            Modulator modulator = new Modulator();

            SpeexDecoder speexDecoder = new SpeexDecoder();
            speexDecoder.init(0, 8000, 1, true);

            AES aes = new AES("abc");
            // byte[] encodedData = "Вставай страна народная".getBytes("UTF-8");
            boolean first = true;
                byte[] encodedData = null;
            while (true) {
                line.read(inBuffer, 0, inBuffer.length);
                speexEncoder.processData(inBuffer, 0, inBuffer.length);
                if (first) {
                    encodedData = new byte[speexEncoder.getProcessedDataByteSize()*2];
                    speexEncoder.getProcessedData(encodedData, 0);
                } else {
                    speexEncoder.getProcessedData(encodedData, encodedData.length/2);
                }




                //   speexDecoder.processData(encodedData,0,encodedData.length);
                //    byte[] decoded = new byte[speexDecoder.getProcessedDataByteSize()];
                //    speexDecoder.getProcessedData(decoded,0);
                //
                //   line2.write(decoded,0,decoded.length);

                if(!first){
                    byte[] ctypted = aes.encrypt(encodedData);
                    byte[] modulated = modulator.modulate(ctypted);
                    AudioChannel.getOutLine().write(modulated, 0, modulated.length);
                    first = true;
                }else{
                    first = false;
                }



            }
        } catch (Throwable t) {
            System.out.println(t);
        }
    }
}