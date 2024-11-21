clear 
clc
close all


% Σταθερό βήμα Ts:
Ts = 0.002;

% Θεωρούμε τον συνολικό χρόνο μέχρι τα 10 sec, ώστε να διαστασιοποιήσομε εξ' αρχής τους πίνακες:
t = 0:Ts:10;

% Ακτίνα μπάλας
r = 0.02;

% Σταθερές Κp και Kv της εισόδου ελέγχου:
Kp = 200;
Kv = 100;

% Επιθυμητοί προσανατολισμοί:
g_be = [cos(pi) 0 sin(pi) 0 ; 0 1 0 0 ; -sin(pi) 0 cos(pi) 0.06 ; 0 0 0 1];
g_oc = [1 0 0 0.4 ; 0 1 0 0 ; 0 0 1 0.2 ; 0 0 0 1];



% Δημιουργία του ρομπότ μέσω της κλήσης του αρχείου mdl_ur10e():
robot = mdl_ur10e();
q1 = [-0.140 -1.556 -1.359 1.425 -1.053 -1.732];




% Κλήση του αρχείου Wspace.p
Wsp = Wspace();


% Αρχικοποίηση πινάκων:
u = zeros(1,6);
Q = zeros(length(t),6);
Q(1,:) = q1;
Q_dot = zeros(length(t),6); 
Q_dot_dot = zeros(length(t),6);
thita_cb = zeros(1, length(t));
thita_cb(1) = -0.3621; %(rad) --> το γνωρίζουμε έπειτα από υπολογισμό
P_be_real = zeros(length(t),3);
k_be_real = zeros(length(t),3);
thita_be_real = zeros(length(t),1);



% Ανάθεση τιμής του ανεκτού σφάλματος:
tol_error = 0.002;



k = 1;
flag_counter = 0;
while(flag_counter < 1/Ts)
    % Κλήση της Wspace:
    [p_cb, v_cb, w_cb] = Wsp.sim_ball(Ts);


    % Υπολογισμός του διανύσματος thita_cb για κάθε χρονική στιγμή:
    if (k<length(t))
        thita_cb(k+1) = thita_cb(k) + w_cb(1)*Ts;
    end


    % Δημιουργία του μετασχηματισμού g_cb κάθε χρονική στιγμή:
    % R_cb = Rot(X, thita_cb)
    g_cb = [1 0 0 p_cb(1) ; 0 cos(thita_cb(k)) -sin(thita_cb(k)) p_cb(2) ; 0 sin(thita_cb(k)) cos(thita_cb(k)) p_cb(3) ; 0 0 0 1];

    % Επιθυμητός προσανατολισμός του άκρου e ως προς την αρχή Ο:
    g_oe_d = g_oc * g_cb * g_be;
    
    % Επιθυμητή θέση του άκρου e ως προς το Ο:
    p_oe_d = [g_oe_d(1,4) g_oe_d(2,4) g_oe_d(3,4)];
    
    % Πραγματικός προσανατολισμός του άκρου e ως προς την αρχή Ο:
    g_oe = robot.fkine(Q(k,:));
    g_oe = g_oe.T;
    p_oe = [g_oe(1,4), g_oe(2,4), g_oe(3,4)];


    % ===== Για την εύρεση του πραγματικού προσανατολισμού g_be =====
    g_be_real = inv(g_cb) * inv(g_oc) * g_oe;
    % Υπολογισμός της θέσης p_be_real:
    P_be_real(k,:) = [g_be_real(1,4) g_be_real(2,4) g_be_real(3,4)];
    % Υπολογισμός του πίνακα R_be_real:
    R_be_real = [g_be_real(1,1) g_be_real(1,2) g_be_real(1,3) ; g_be_real(2,1) g_be_real(2,2) g_be_real(2,3) ; g_be_real(3,1) g_be_real(3,2) g_be_real(3,3)];
    % Υπολογισμός του ισοδύναμου άξονα-γωνίας του R_be_real:
    Trace_be_real = R_be_real(1,1) +  R_be_real(2,2) +  R_be_real(3,3);
    thita_be_real(k) = acos((Trace_be_real-1)/2);
    k_be_real(k,:) = (1/(2*sin(thita_be_real(k))))*[R_be_real(3,2)-R_be_real(2,3) ; R_be_real(1,3)-R_be_real(3,1) ; R_be_real(2,1)-R_be_real(1,2)];
    % Κανονικοποίηση του διανύσματος k_be_real ώστε να έχει μέτρο ίσο με 1:
    if(norm(k_be_real(k,:))>1)
        k_be_real(k,:) = k_be_real(k,:)/norm(k_be_real(k,:));
    end


    % Υπολογισμός Ιακωβιανής του βραχίονα
    J = robot.jacob0(Q(k,:));

    
    R_oe = [g_oe(1,1) g_oe(1,2) g_oe(1,3) ; g_oe(2,1) g_oe(2,2) g_oe(2,3) ; g_oe(3,1) g_oe(3,2) g_oe(3,3)];
    R_oe_d = [g_oe_d(1,1) g_oe_d(1,2) g_oe_d(1,3) ; g_oe_d(2,1) g_oe_d(2,2) g_oe_d(2,3) ; g_oe_d(3,1) g_oe_d(3,2) g_oe_d(3,3)];
    % Σφάλμα προσανατολισμού R_e:
    R_e = R_oe * R_oe_d';
    Trace = R_e(1,1) + R_e(2,2) + R_e(3,3);
    thita_e = acos((Trace-1)/2);
    ke = (1/(2*sin(thita_e)))*[R_e(3,2)-R_e(2,3) ; R_e(1,3)-R_e(3,1) ; R_e(2,1)-R_e(1,2)];
    % Κανονικοποίηση του διανύσματος ke ώστε να έχει μέτρο ίσο με 1:
    if(norm(ke)>1)
        ke = ke/norm(ke);
    end

    % Σφάλμα προσανατολισμού e_l:
    e_l = thita_e * ke;


    % Υπολογισμός της εισόδου ελέγχου u του συστήματος:
    u(1) = v_cb(1) - Kp * (p_oe(1) - p_oe_d(1));
    u(2) = v_cb(2) - Kp * (p_oe(2)-p_oe_d(2));
    u(3) = v_cb(3) - Kp * (p_oe(3)-p_oe_d(3));
    u(4) = w_cb(1) - Kv * e_l(1);
    u(5) = w_cb(2) - Kv * e_l(2);
    u(6) = w_cb(3) - Kv * e_l(3);


    % Υπολογισμός πίνακα ταχυτήτων Q_dot(k,:)
    Q_dot(k,:) = J\u';
    
    % Πέρασμα από συνάρτηση κορεσμού για τα όρια ταχύτητας:
    Q_dot(k,1) = sat(Q_dot(k,1), -(2/3)*pi, (2/3)*pi);
    Q_dot(k,2) = sat(Q_dot(k,2), -(2/3)*pi, (2/3)*pi);
    Q_dot(k,3) = sat(Q_dot(k,3), -pi, pi);
    Q_dot(k,4) = sat(Q_dot(k,4), -pi, pi);
    Q_dot(k,5) = sat(Q_dot(k,5), -pi, pi);
    Q_dot(k,6) = sat(Q_dot(k,6), -pi, pi);
    

    % Υπολογισμός του πίνακα επιταχύνσεων Q_dot_dot(k,:)
    if(k>1)
        Q_dot_dot(k-1,:) = (Q_dot(k,:)-Q_dot(k-1,:))/Ts;
        % Συνάρτηση κορεσμού για τα όρια επιτάχυνσης:
        Q_dot_dot(k-1,1) = sat(Q_dot_dot(k-1,1), -250, 250);
        Q_dot_dot(k-1,2) = sat(Q_dot_dot(k-1,2), -250, 250);
        Q_dot_dot(k-1,3) = sat(Q_dot_dot(k-1,3), -250, 250);
        Q_dot_dot(k-1,4) = sat(Q_dot_dot(k-1,4), -250, 250);
        Q_dot_dot(k-1,5) = sat(Q_dot_dot(k-1,5), -250, 250);
        Q_dot_dot(k-1,6) = sat(Q_dot_dot(k-1,6), -250, 250);

        % Εκ νέου υπολογισμός του Q_dot(k,:) ώστε να πληρούνται όλα τα όρια:
        Q_dot(k,:) = Q_dot(k-1,:) + Q_dot_dot(k-1,:)*Ts;
    end

    
    % Υπολογισμός του πίνακα θέσεων Q(k+1)
    if(k<length(t))
        Q(k+1,:) = Q(k,:) + Q_dot(k,:) * Ts;
    end

    

    % Σφάλμα θέσης:
    e_p = p_oe - p_oe_d;

    % Έλεγχος εαν σε αυτό το βήμα τα σφάλματα θέσης και προσανατολισμού είναι μικρότερα από
    % το ανεκτό σφάλμα tol_error:
    if (abs(e_l(1)) < tol_error && abs(e_l(2)) < tol_error && abs(e_l(3)) < tol_error ...
            && abs(e_p(1)) < tol_error && abs(e_p(2)) < tol_error && abs(e_p(3)) < tol_error)
        flag_counter = flag_counter + 1;
    else
        flag_counter = 0;
    end

    
    k = k+1;
end






% Περιορισμός του χρόνου t μέχρι τη χρονική στιγμή επίτευξης του στόχου:
t = t(1:k);
% Περιορισμός των πινάκων μέχρι τη χρονική στιγμή επίτευξης του στόχου:
Q = Q(1:k, 1:6);
Q_dot = Q_dot(1:k-1, 1:6);
Q_dot_dot = Q_dot_dot(1:k-1, 1:6);
P_be_real = P_be_real(1:k-1, 1:3);
k_be_real = k_be_real(1:k-1, 1:3);
thita_be_real = thita_be_real(1:k-1);





% Αποθήκευση γραφημάτων:
% ========== Θέσεις αρθρώσεων ===========
figure(1);
counter = 1;
for i=1:1:6
    subplot(2,3,i)
    plot(t,Q(:,i));
    xlim([0, t(k)])
    xlabel("t (sec)")
    ylabel("q_" +i+ " (rad)")
    title("Απόκριση q_" +i)
end

% ========== Ταχύτητες αρθρώσεων ===========
figure(2);
for i=1:1:6
    subplot(2,3,i)
    plot(t(1:k-1),Q_dot(:,i));
    xlim([0, t(k-1)])
    xlabel("t (sec)")
    ylabel("qdot_" +i+ " (rad/sec)")
    title("Απόκριση qdot_" +i)
end

% ========== Επιταχύνσεις αρθρώσεων ===========
figure(3);
for i=1:1:6
    subplot(2,3,i)
    plot(t(1:k-1),Q_dot_dot(:,i));
    xlim([0, t(k-1)])
    xlabel("t (sec)")
    ylabel("qdotdot_" +i+ " (rad/sec^2)")
    title("Απόκριση qdotdot_" +i)
end

% ======= Πραγματική θέση p_be_real ========
figure(4);
subplot(3,1,1)
plot(t(1:k-1),P_be_real(:,1))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("pbe(x)  (meter)")
title("Πραγματική θέση pbe real_x")

subplot(3,1,2)
plot(t(1:k-1),P_be_real(:,2))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("pbe(y)  (meter)")
title("Πραγματική θέση pbe real_y")

subplot(3,1,3)
plot(t(1:k-1),P_be_real(:,3))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("pbe(z)  (meter)")
title("Πραγματική θέση pbe real_z")

% ===== R_be_real ισοδύναμος άξονας k_be_real =====
figure(5);
subplot(3,1,1)
plot(t(1:k-1), k_be_real(:,1))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("k be_x")
title("k be_x")

subplot(3,1,2)
plot(t(1:k-1), k_be_real(:,2))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("k be_y")
title("k be_y")

subplot(3,1,3)
plot(t(1:k-1), k_be_real(:,3))
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("k be_z")
title("k be_z")

% ===== R_be_real ισοδύναμη γωνία thita_be_real =====
figure(6)
plot(t(1:k-1), thita_be_real)
xlim([0, t(k-1)])
xlabel("t (sec)")
ylabel("thita be real  (rad)")
title("thita be real")


% ==== Οπτικοποίηση της παρακολούθησης της μπάλας από το ρομπότ: ====
figure(7)
view = [90 0];
Wsp.visualize(robot,Q',view);







% Συνάρτηση κορεσμού:
function f = sat(input, lower_limit, upper_limit)
    f = min(upper_limit, max(lower_limit, input));
end


