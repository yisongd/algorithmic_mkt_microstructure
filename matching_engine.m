%%% parameters %%%
buy = 1;
sell = 2;
t_min = 1;
t_max = 6323;
burn_in_period = 3322;
price_max = 20;
price_min = 1;
volum_max = 100;
volum_min = 1;
num_bgt = 10;
order_id = 0;

%%% data structure %%%
current_entry = zeros(1,7);

live_buy_orders_list=zeros(t_max,7);
live_sell_orders_list=ones(t_max,7)*1000;

transaction_price_volume_stor_mat=zeros(1,7);

best_bid_history = zeros(t_max,2);
best_ask_history = zeros(t_max,2);
bid_ask_stor_mat=zeros(t_max,2);
bid_ask_depth_stor_mat = zeros(t_max,2);
bid_ask_spread = zeros(t_max,2);
LOB=zeros(t_max,2);

z = 3;
robot_z_cash_position=zeros(t_max,3);
robot_z_trading_profit=zeros(t_max,3);
robot_z_inventory = zeros(t_max,3);
robot_z_inventory_position=zeros(t_max,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;
tic
while t <= t_max
    
    
    %%%%%%% generate random limit order %%%%%%
    order_id = order_id + 1;
    buy_sell_indi = randi([0,1],1,1);
    if buy_sell_indi == 0
        buy_sell_indi = -1;
    end
    
    current_entry = [randi([1,num_bgt], 1,1), buy_sell_indi, randi([price_min,price_max],1,1), randi([volum_min,volum_max],1,1), order_id, t, 1];
    
    if buy_sell_indi == 1
        live_buy_orders_list(t,:) = current_entry;
        live_buy_orders_list = sortrows(live_buy_orders_list,[-3,6]);
    end
    if buy_sell_indi == -1
        live_sell_orders_list(t,:) = current_entry;
        live_sell_orders_list = sortrows(live_sell_orders_list,[3,6]); 
    end
    %%%%%%%%% generate random limit order finish %%%%%%%%%%
    
    
    
    %%%%%%%% matching engine algo %%%%%%%% 
    if buy_sell_indi == 1
        aggressor_sign = buy_sell_indi;
        aggressor_account_id = live_buy_orders_list(1,1);
        %if the buying entry can be traded
        while live_buy_orders_list(1,3)>=live_sell_orders_list(1,3) && live_buy_orders_list(1,4)>=live_sell_orders_list(1,4) 
            if live_sell_orders_list(1,1)==live_buy_orders_list(1,1)      %prevent self-transaction
                live_sell_orders_list(1,:) = ones(1,7)*1000;
                live_sell_orders_list = sortrows(live_sell_orders_list, [3,6]);
                continue
            end
            passor_account_id = live_sell_orders_list(1,1);
            passor_order_id = live_sell_orders_list(1,5);
            transac_quant = live_sell_orders_list(1,4);
            transac_price = live_sell_orders_list(1,3);
            transaction_price_volume_stor_mat = [transaction_price_volume_stor_mat;[t, aggressor_sign, transac_price, transac_quant, passor_order_id, passor_account_id,aggressor_account_id]];  
            live_buy_orders_list(1,4) = live_buy_orders_list(1,4) - live_sell_orders_list(1,4);
            if live_buy_orders_list(1,4)==0
                live_sell_orders_list(1,:) = ones(1,7)*1000;
                live_sell_orders_list = sortrows(live_sell_orders_list,[3,6]);
                live_buy_orders_list(1,:) = zeros(1,7);
                live_buy_orders_list = sortrows(live_buy_orders_list,[-3,6]);
                break
            end
            live_sell_orders_list(1,:) = ones(1,7)*1000;
            live_sell_orders_list = sortrows(live_sell_orders_list,[3,6]);
        end
        while live_buy_orders_list(1,4)<live_sell_orders_list(1,4) && live_buy_orders_list(1,3)>=live_sell_orders_list(1,3)
            if live_sell_orders_list(1,1)==live_buy_orders_list(1,1)
                live_sell_orders_list(1,:) = ones(1,7)*1000;
                live_sell_orders_list = sortrows(live_sell_orders_list, [3,6]);
                continue
            end
            passor_account_id = live_sell_orders_list(1,1);
            passor_order_id = live_sell_orders_list(1,5);
            transac_price = live_sell_orders_list(1,3);
            transac_quant = live_buy_orders_list(1,4);
            transaction_price_volume_stor_mat = [transaction_price_volume_stor_mat;[t, aggressor_sign, transac_price, transac_quant, passor_order_id, passor_account_id,aggressor_account_id]];  
            live_sell_orders_list(1,4) = live_sell_orders_list(1,4) - live_buy_orders_list(1,4);
            live_buy_orders_list(1,:) = zeros(1,7);
            live_buy_orders_list = sortrows(live_buy_orders_list,[-3,6]);
        end
    end
    
    if buy_sell_indi == -1
        aggressor_sign = buy_sell_indi;
        aggressor_account_id = live_sell_orders_list(1,1);
        %if the selling entry can be traded
        while live_sell_orders_list(1,3)<=live_buy_orders_list(1,3) && live_sell_orders_list(1,4)>=live_buy_orders_list(1,4)
            if live_sell_orders_list(1,1)==live_buy_orders_list(1,1) %prevent self-transaction
                live_buy_orders_list(1,:) = zeros(1,7);
                live_buy_orders_list = sortrows(live_buy_orders_list, [-3,6]);
                continue
            end
            passor_account_id = live_buy_orders_list(1,1);
            passor_order_id = live_buy_orders_list(1,5);
            transac_quant = live_buy_orders_list(1,4);
            transac_price = live_buy_orders_list(1,3);
            transaction_price_volume_stor_mat = [transaction_price_volume_stor_mat; [t, aggressor_sign, transac_price, transac_quant, passor_order_id, passor_account_id,aggressor_account_id]];  
            live_sell_orders_list(1,4) = live_sell_orders_list(1,4) - live_buy_orders_list(1,4);
            if live_sell_orders_list(1,4)==0
                live_sell_orders_list(1,:) = ones(1,7)*1000;
                live_sell_orders_list = sortrows(live_sell_orders_list,[3,6]);
                live_buy_orders_list(1,:) = zeros(1,7);
                live_buy_orders_list = sortrows(live_buy_orders_list,[-3,6]);
                break
            end
            live_buy_orders_list(1,:) = zeros(1,7);
            live_buy_orders_list = sortrows(live_buy_orders_list,[-3,6]);
        end
        while live_sell_orders_list(1,4)<live_buy_orders_list(1,4) && live_buy_orders_list(1,3)>=live_sell_orders_list(1,3)
            if live_sell_orders_list(1,1)==live_buy_orders_list(1,1)
                live_buy_orders_list(1,:) = zeros(1,7);
                live_buy_orders_list = sortrows(live_buy_orders_list, [-3,6]);
                continue
            end
            passor_account_id = live_buy_orders_list(1,1);
            passor_order_id = live_buy_orders_list(1,5);
            transac_price = live_buy_orders_list(1,3);
            transac_quant = live_sell_orders_list(1,4);
            transaction_price_volume_stor_mat = [transaction_price_volume_stor_mat;[t, aggressor_sign, transac_price, transac_quant, passor_order_id, passor_account_id,aggressor_account_id]];
            live_buy_orders_list(1,4) = live_buy_orders_list(1,4) - live_sell_orders_list(1,4);
            live_sell_orders_list(1,:) = ones(1,7)*1000;
            live_sell_orders_list = sortrows(live_sell_orders_list,[3,6]);
        end
    end
    %%%%%%%%%%%%%%%%%%matching engine algo finish %%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%% constructing orderbook depth %%%%%
    
    % current best bid and best ask
    best_bid = live_buy_orders_list (1,3);
    best_ask = live_sell_orders_list (1,3);
    
    % current best bid depth and best ask depth
    best_bid_depth = 0;
    best_ask_depth = 0;
    i = 1;
    while i<t_max
        if live_buy_orders_list(i,3)==live_buy_orders_list(1,3) && live_buy_orders_list(1,3)~=0
            best_bid_depth = best_bid_depth + live_buy_orders_list(1,4);
        end
        if live_sell_orders_list(i,3)==live_sell_orders_list(1,3) && live_sell_orders_list(1,3)~=1000
            best_ask_depth = best_ask_depth + live_sell_orders_list(1,4);
        end
        i = i + 1;
    end
    
    % restore all best bid depth & best ask depth
    bid_ask_depth_stor_mat(t,:) = [best_bid_depth,best_ask_depth];
    
    
    % restore all bid-ask spread
    if live_buy_orders_list(1,3)~=0 && live_sell_orders_list(1,3)~=1000
        bid_ask_stor_mat(t,:) = [live_buy_orders_list(1,3),live_sell_orders_list(1,3)];
    end
    bid_ask_spread(t,:) = [t,bid_ask_stor_mat(t,2)-bid_ask_stor_mat(t,1)];
    
    %update live buy depth
    price_quant_of_live_buy = live_buy_orders_list(:,[3,4]);
    [a,~,c]=unique(price_quant_of_live_buy(:,1));
    buy_lob_depth = [a, accumarray(c,price_quant_of_live_buy(:,2))];
    buy_lob_depth(~any( buy_lob_depth, 2), : ) = [ ]; %[price, total quant]

    %update live sell depth
    price_quant_of_live_sell = live_sell_orders_list(:,[3,4]);
    [a,~,c]=unique(price_quant_of_live_sell(:,1));
    sell_lob_depth = [a, accumarray(c,price_quant_of_live_sell(:,2))];
    indicator_1000 = sell_lob_depth(:,1)==1000;
    if numel(indicator_1000>1)
        sell_lob_depth(indicator_1000(1:end),:)=[]; %[price, total quant]
    end

    %combine live sell depth & live buy depth, covenient for plotting
    LOB = [buy_lob_depth;sell_lob_depth];    
    
    %%%%%%%%%%%% constructing orderbook depth finish %%%%%%%%%%%%%%%
    t=t+1;
end

% delete all empty rows in the out put matrix
live_buy_orders_list(~any(live_buy_orders_list,2),:)=[]; % in live_buy_orders
indi1000 = live_sell_orders_list(:,1)==1000; % in live_sell_orders
if numel(indi1000>1)
    live_sell_orders_list(indi1000(1:end),:)=[];
end
toc

LOB(~any(LOB,2),:)=[]; % in combining depth matrix

%%%%% robot_z_cash_inventory_profit %%%%%
i_1 = 1;
current_cash_position = 0;
current_inventory = 0;
[m,n] = size(transaction_price_volume_stor_mat);
for i_1 = 1:m
    % if z is aggressor
    if transaction_price_volume_stor_mat(i_1,7) == z
        current_cash_position = current_cash_position - transaction_price_volume_stor_mat(i_1,2)* transaction_price_volume_stor_mat(i_1,3)*transaction_price_volume_stor_mat(i_1,4);
        current_inventory = current_inventory + transaction_price_volume_stor_mat(i_1,2)*transaction_price_volume_stor_mat(i_1,4);
    end
    % if z is passor
    if transaction_price_volume_stor_mat(i_1,6) == z
        current_cash_position = current_cash_position + transaction_price_volume_stor_mat(i_1,2)*transaction_price_volume_stor_mat(i_1,3)*transaction_price_volume_stor_mat(i_1,4);
        current_inventory = current_inventory - transaction_price_volume_stor_mat(i_1,2)*transaction_price_volume_stor_mat(i_1,4);
    end
    current_inventory_value = current_inventory*transaction_price_volume_stor_mat(i_1,3);
    trading_profit = current_cash_position + current_inventory_value;
    robot_z_cash_position(i_1,:) = [z,i_1,current_cash_position];
    robot_z_inventory_position(i_1,:) = [z,i_1,current_inventory_value];
    robot_z_inventory(i_1,:) = [z,i_1,current_inventory];
    robot_z_trading_profit(i_1,:) = [z,i_1,trading_profit];
end
robot_z_cash_position(  ~any( robot_z_cash_position, 2), : ) = [ ];
robot_z_inventory_position(  ~any( robot_z_inventory_position, 2), : ) = [ ];
robot_z_trading_profit(  ~any( robot_z_trading_profit, 2), : ) = [ ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% analysis %%%%
total_number_of_trades = size(transaction_price_volume_stor_mat);

subplot(3,2,1);
plot(bid_ask_spread(:,1),bid_ask_spread(:,2));
ylim([0 price_max-1]);
xlim([0 t_max]);

subplot(3,2,2);

subplot(3,2,3);
bar(LOB(:,1),LOB(:,2));

subplot(3,2,4);
plot(bid_ask_stor_mat(:,1),'r');
hold on
plot(bid_ask_stor_mat(:,2),'y');
ylim([ 0 price_max]);

subplot(3,2,5);
bar(bid_ask_depth_stor_mat(:,1),'r');
hold on
bar(bid_ask_depth_stor_mat(:,2),'y');

subplot(3,2,6);
plot(robot_z_cash_position(:,2),robot_z_cash_position(:,3));
hold on
plot(robot_z_inventory_position(:,2),robot_z_inventory_position(:,3),'y')
hold on
plot(robot_z_trading_profit(:,2),robot_z_trading_profit(:,3),'r')
% subplot(4,2,6);
% plot(robot_z_cash_position(:,2),robot_z_cash_position(:,3));
% 
% subplot(4,2,7);
% plot(robot_z_inventory_position(:,2),robot_z_inventory_position(:,3));
% 
% subplot(4,2,8);
% plot(robot_z_trading_profit(:,2),robot_z_trading_profit(:,3));