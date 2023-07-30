import Icon from './svg/immune_system.svg';

interface Icon {
  size: number;
}

const ImmuneIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="immune system icon" />
    </div>
  );
};

export default ImmuneIcon;
