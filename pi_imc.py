import numpy as np
import matplotlib.pyplot as plt
from control import tf, feedback, margin, bode_plot, frequency_response, step_response, bandwidth, pade

class PIIMC:
    def __init__(self, gain, tau, theta, padeDeg=2, *args):
        self.gain = gain
        self.tau = tau
        self.theta = theta
        s = tf('s')
        delayApprox = pade(theta, padeDeg)
        sysd = tf(delayApprox[0], delayApprox[1])
        self.P = sysd*gain / (tau*s+1)
        #self.P = np.exp(-theta * s) * gain / (tau * s + 1)  # no explict time delay support?        
        print(args[0])
        if len(args) == 1:
            self.calcpigains(args[0])
        elif len(args) == 2:
            self.setpigains(args[0], args[1])
        elif len(args) < 1 or len(args) > 2:
            raise ValueError('Constructor takes 3, 4, or 5 arguments')

    def calcpigains(self, fac):
        '''Calculate PI gains based on IMC tuning factor'''
        self.Tc = max(fac * 0.80 * self.theta, fac * 0.1 * self.tau)
        Kc = (1 / self.gain) * (self.tau / (self.theta + self.Tc))
        Ti = self.tau
        self.Kp = Kc
        self.Ki = Kc / Ti
        s = tf('s')
        self.C = (self.Kp * s + self.Ki) / s
        self.So = 1 / (self.C * self.P + 1)
        self.Si = self.P / (self.C * self.P + 1)
        gm_, pm_, wcg, wcp = margin(self.C * self.P)
        self.Gm = 20 * np.log10(gm_)
        self.Pm = pm_
        self.wcg = wcg
        self.wcp = wcp
        self.BWo = bandwidth(feedback(self.C * self.P, 1))
        self.To = 1 / (self.BWo / (2 * np.pi))

    def setpigains(self, Kp, Ki):
        '''Calculate trans funcs. based on user chosen PI gains'''
        self.Kp = Kp
        self.Ki = Ki
        s = tf('s')
        self.C = (self.Kp * s + self.Ki) / s
        self.So = 1 / (self.C * self.P + 1)
        self.Si = self.P / (self.C * self.P + 1)
        gm_, pm_, wcg, wcp = margin(self.C * self.P)
        self.Gm = 20 * np.log10(gm_)
        self.Pm = pm_
        self.wcg = wcg
        self.wcp = wcp
        self.BWo = bandwidth(feedback(self.C * self.P, 1))
        self.To = 1 / (self.BWo / (2 * np.pi))

    def margins(self):
        '''Gain and phase margins for the loop function: CP)'''
        plt.figure()        
        figTitle = 'Loop Gain & Phase Margins'
        bode_plot(self.C * self.P, display_margins=True, title=figTitle)
        plt.show()

    def step(self):
        '''Step response plot for closed-loop trans func: CP/(1+CP)'''
        plt.figure()
        t, y = step_response(feedback(self.C * self.P, 1))
        plt.plot(t, y)
        plt.grid(True)
        plt.title('Closed-Loop Step Response')
        plt.xlabel('Time')
        plt.ylabel('Output')
        plt.show()

    def closedloopbode(self):
        '''Bode plot for closed-loop trans func: CP/(1+CP)'''
        plt.figure()        
        mag, phase, omega = frequency_response(feedback(self.C * self.P, 1))
        plt.subplot(211)
        plt.semilogx(omega, 20 * np.log10(mag))
        plt.title('Closed Loop Bode Diagram')
        plt.ylabel('Gain [dB]')
        plt.grid(True)
        plt.subplot(212)
        plt.semilogx(omega, phase * 180 / np.pi)
        plt.ylabel('Phase [deg]')
        plt.xlabel('Freq. [rad/s]')
        plt.grid(True)
        plt.show()

    def Sobode(self):
        '''Bode plot for output disturbance trans func: 1/(1+CP)'''
        plt.figure()
        mag, phase, omega = frequency_response(1 / (1 + self.C * self.P))
        plt.subplot(211)
        plt.semilogx(omega, 20 * np.log10(mag))
        plt.title('Output Dist. Bode Diagram')
        plt.ylabel('Gain [dB]')
        plt.grid(True)
        plt.subplot(212)
        plt.semilogx(omega, phase * 180 / np.pi)
        plt.ylabel('Phase [deg]')
        plt.xlabel('Freq. [rad/s]')
        plt.grid(True)
        plt.show()

    def Sibode(self):
        '''Bode plot for input disturbance trans func: P/(1+CP)'''
        plt.figure()
        mag, phase, omega = frequency_response(self.P / (1 + self.C * self.P))
        plt.subplot(211)
        plt.semilogx(omega, 20 * np.log10(mag))
        plt.title('Input Dist. Bode Diagram')
        plt.ylabel('Gain [dB]')
        plt.grid(True)
        plt.subplot(212)
        plt.semilogx(omega, phase * 180 / np.pi)
        plt.ylabel('Phase [deg]')
        plt.xlabel('Freq. [rad/s]')
        plt.grid(True)
        plt.show()

# Example usage:
# pi_controller = PIIMC(gain=1.0, tau=1.0, theta=0.5, pade =1, fac=1.5)
# pi_controller.step()