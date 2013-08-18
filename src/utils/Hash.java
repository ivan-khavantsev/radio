package utils;

import java.io.UnsupportedEncodingException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class Hash {
    public static String md5(String message) throws UnsupportedEncodingException, NoSuchAlgorithmException {
        String digest = null;

        MessageDigest md = MessageDigest.getInstance("MD5");
        byte[] hash = md.digest(message.getBytes("UTF-8")); //converting byte array to Hexadecimal
        StringBuilder sb = new StringBuilder(2 * hash.length);
        for (byte b : hash) {
            sb.append(String.format("%02x", b & 0xff));
        }
        digest = sb.toString();

        return digest;
    }

    public static byte[] md5(byte[] message) throws NoSuchAlgorithmException {
        MessageDigest md = MessageDigest.getInstance("MD5");
        byte[] hash = md.digest(message);
        return hash;
    }

}
