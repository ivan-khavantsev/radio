package crypt;

import utils.Hash;

import javax.crypto.*;
import javax.crypto.spec.SecretKeySpec;
import java.io.UnsupportedEncodingException;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;
import java.security.NoSuchProviderException;

public class AES {
    private Cipher cipher = null;
    SecretKey key = null;
    private static String ALGORYTHM = "AES/ECB/PKCS5Padding";
    private static String PROVIDER = "SunJCE";

    public AES(String password) throws NoSuchAlgorithmException, NoSuchProviderException, NoSuchPaddingException, UnsupportedEncodingException {
        byte[] passwordHash = Hash.md5(password.getBytes("UTF-8"));
        this.key = new SecretKeySpec(passwordHash, "AES");
        this.cipher = Cipher.getInstance(ALGORYTHM, PROVIDER);
    }

    public byte[] encrypt(byte[] data) throws
            InvalidKeyException, IllegalBlockSizeException, BadPaddingException {

        cipher.init(Cipher.ENCRYPT_MODE, key);
        byte[] encrypted = cipher.doFinal(data);
        return encrypted;
    }

    public byte[] decrypt(byte[] encypted) throws
            InvalidKeyException, IllegalBlockSizeException, BadPaddingException {

        cipher.init(Cipher.DECRYPT_MODE, key);
        byte[] decrypted = cipher.doFinal(encypted);
        return decrypted;
    }

    public static void main(String[] args) {
        String text = "Some text";
        String password = "abcd";
        try {
            System.out.println("text is " + text);
            System.out.println("password is " + password);
            AES aes = new AES(password);
            byte[] cripted = aes.encrypt(text.getBytes("UTF-8"));
            System.out.println("crypt is " + new String(cripted));
            String decripted = new String(aes.decrypt(cripted));
            System.out.println("decrypted is " + decripted);
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }
}