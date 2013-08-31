package audio;

import javax.sound.sampled.*;

public class AudioChannel {
    private static AudioFormat channelAudioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, 44100, 16, 1, 2, 44100, false);
    private static SourceDataLine out = null;
    private static TargetDataLine in = null;

    static {
        try {
            out = AudioSystem.getSourceDataLine(channelAudioFormat);
            out.open();
            out.start();

            in = (TargetDataLine) AudioSystem.getLine(new DataLine.Info(TargetDataLine.class, channelAudioFormat));
            in.open();
            in.start();

        } catch (Throwable t) {
            System.out.println(t);
        }

    }

    public static TargetDataLine getInLine() {
        return in;
    }

    public static SourceDataLine getOutLine() {
        return out;
    }

}
